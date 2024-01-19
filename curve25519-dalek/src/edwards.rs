// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2020 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - isis agora lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>

//! Group operations for Curve25519, in Edwards form.
//!
//! ## Encoding and Decoding
//!
//! Encoding is done by converting to and from a `CompressedEdwardsY`
//! struct, which is a typed wrapper around `[u8; 32]`.
//!
//! ## Equality Testing
//!
//! The `EdwardsPoint` struct implements the [`subtle::ConstantTimeEq`]
//! trait for constant-time equality checking, and the Rust `Eq` trait
//! for variable-time equality checking.
//!
//! ## Cofactor-related functions
//!
//! The order of the group of points on the curve \\(\mathcal E\\)
//! is \\(|\mathcal E| = 8\ell \\), so its structure is \\( \mathcal
//! E = \mathcal E\[8\] \times \mathcal E[\ell]\\).  The torsion
//! subgroup \\( \mathcal E\[8\] \\) consists of eight points of small
//! order.  Technically, all of \\(\mathcal E\\) is torsion, but we
//! use the word only to refer to the small \\(\mathcal E\[8\]\\) part, not
//! the large prime-order \\(\mathcal E[\ell]\\) part.
//!
//! To test if a point is in \\( \mathcal E\[8\] \\), use
//! [`EdwardsPoint::is_small_order`].
//!
//! To test if a point is in \\( \mathcal E[\ell] \\), use
//! [`EdwardsPoint::is_torsion_free`].
//!
//! To multiply by the cofactor, use [`EdwardsPoint::mul_by_cofactor`].
//!
//! To avoid dealing with cofactors entirely, consider using Ristretto.
//!
//! ## Scalars
//!
//! Scalars are represented by the [`Scalar`] struct. To construct a scalar, see
//! [`Scalar::from_canonical_bytes`] or [`Scalar::from_bytes_mod_order_wide`].
//!
//! ## Scalar Multiplication
//!
//! Scalar multiplication on Edwards points is provided by:
//!
//! * the `*` operator between a `Scalar` and a `EdwardsPoint`, which
//! performs constant-time variable-base scalar multiplication;
//!
//! * the `*` operator between a `Scalar` and a
//! `EdwardsBasepointTable`, which performs constant-time fixed-base
//! scalar multiplication;
//!
//! * an implementation of the
//! [`MultiscalarMul`](../traits/trait.MultiscalarMul.html) trait for
//! constant-time variable-base multiscalar multiplication;
//!
//! * an implementation of the
//! [`VartimeMultiscalarMul`](../traits/trait.VartimeMultiscalarMul.html)
//! trait for variable-time variable-base multiscalar multiplication;
//!
//! ## Implementation
//!
//! The Edwards arithmetic is implemented using the “extended twisted
//! coordinates” of Hisil, Wong, Carter, and Dawson, and the
//! corresponding complete formulas.  For more details,
//! see the [`curve_models` submodule][curve_models]
//! of the internal documentation.
//!
//! ## Validity Checking
//!
//! There is no function for checking whether a point is valid.
//! Instead, the `EdwardsPoint` struct is guaranteed to hold a valid
//! point on the curve.
//!
//! We use the Rust type system to make invalid points
//! unrepresentable: `EdwardsPoint` objects can only be created via
//! successful decompression of a compressed point, or else by
//! operations on other (valid) `EdwardsPoint`s.
//!
//! [curve_models]: https://docs.rs/curve25519-dalek/latest/curve25519-dalek/backend/serial/curve_models/index.html

// We allow non snake_case names because coordinates in projective space are
// traditionally denoted by the capitalisation of their respective
// counterparts in affine space.  Yeah, you heard me, rustc, I'm gonna have my
// affine and projective cakes and eat both of them too.
#![allow(non_snake_case)]

use core::array::TryFromSliceError;
use core::borrow::Borrow;
use core::fmt::Debug;
use core::iter::Iterator;
use core::iter::Sum;
use core::ops::{Add, Neg, Sub};
use core::ops::{AddAssign, SubAssign};
use core::ops::{Mul, MulAssign};

use cfg_if::cfg_if;

#[cfg(feature = "digest")]
use digest::{generic_array::typenum::U64, Digest};

#[cfg(feature = "group")]
use {
    group::{cofactor::CofactorGroup, prime::PrimeGroup, GroupEncoding},
    rand_core::RngCore,
    subtle::CtOption,
};

use subtle::Choice;
use subtle::ConditionallyNegatable;
use subtle::ConstantTimeEq;

use crate::constants;

use crate::field::FieldElement;
use crate::scalar::Scalar;

use crate::backend::serial::curve_models::AffineNielsPoint;
use crate::backend::serial::curve_models::CompletedPoint;
use crate::backend::serial::curve_models::ProjectiveNielsPoint;
use crate::backend::serial::curve_models::ProjectivePoint;

#[cfg(feature = "precomputed-tables")]
use crate::window::{
    LookupTableRadix128, LookupTableRadix16, LookupTableRadix256, LookupTableRadix32,
    LookupTableRadix64,
};

#[cfg(feature = "precomputed-tables")]
use crate::traits::BasepointTable;

use crate::traits::{Identity};



// ------------------------------------------------------------------------
// Compressed points
// ------------------------------------------------------------------------

/// In "Edwards y" / "Ed25519" format, the curve point \\((x,y)\\) is
/// determined by the \\(y\\)-coordinate and the sign of \\(x\\).
///
/// The first 255 bits of a `CompressedEdwardsY` represent the
/// \\(y\\)-coordinate.  The high bit of the 32nd byte gives the sign of \\(x\\).
#[derive(Copy, Clone, Eq, PartialEq, Hash)]
pub struct CompressedEdwardsY(pub [u8; 32]);

impl ConstantTimeEq for CompressedEdwardsY {
    fn ct_eq(&self, other: &CompressedEdwardsY) -> Choice {
        self.as_bytes().ct_eq(other.as_bytes())
    }
}

impl Debug for CompressedEdwardsY {
    fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
        write!(f, "CompressedEdwardsY: {:?}", self.as_bytes())
    }
}

impl CompressedEdwardsY {
    /// View this `CompressedEdwardsY` as an array of bytes.
    pub(crate) const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }

    /// Attempt to decompress to an `EdwardsPoint`.
    ///
    /// Returns `None` if the input is not the \\(y\\)-coordinate of a
    /// curve point.
    fn decompress(&self) -> Option<EdwardsPoint> {
        let (is_valid_y_coord, X, Y, Z) = decompress::step_1(self);

        if is_valid_y_coord.into() {
            Some(decompress::step_2(self, X, Y, Z))
        } else {
            None
        }
    }

    pub fn is_curve_point(&self) -> bool {
        self.decompress().is_some()
    }
}

mod decompress {
    use super::*;

    #[rustfmt::skip] // keep alignment of explanatory comments
    pub(super) fn step_1(
        repr: &CompressedEdwardsY,
    ) -> (Choice, FieldElement, FieldElement, FieldElement) {
        let Y = FieldElement::from_bytes(repr.as_bytes());
        let Z = FieldElement::ONE;
        let YY = Y.square();
        let u = &YY - &Z;                            // u =  y²-1
        let v = &(&YY * &constants::EDWARDS_D) + &Z; // v = dy²+1
        let (is_valid_y_coord, X) = FieldElement::sqrt_ratio_i(&u, &v);

        (is_valid_y_coord, X, Y, Z)
    }

    #[rustfmt::skip]
    pub(super) fn step_2(
        repr: &CompressedEdwardsY,
        mut X: FieldElement,
        Y: FieldElement,
        Z: FieldElement,
    ) -> EdwardsPoint {
         // FieldElement::sqrt_ratio_i always returns the nonnegative square root,
         // so we negate according to the supplied sign bit.
        let compressed_sign_bit = Choice::from(repr.as_bytes()[31] >> 7);
        X.conditional_negate(compressed_sign_bit);

        EdwardsPoint {
            X,
            Y,
            Z,
            T: &X * &Y,
        }
    }
}

impl TryFrom<&[u8]> for CompressedEdwardsY {
    type Error = TryFromSliceError;

    fn try_from(slice: &[u8]) -> Result<CompressedEdwardsY, TryFromSliceError> {
        Self::from_slice(slice)
    }
}


// ------------------------------------------------------------------------
// Internal point representations
// ------------------------------------------------------------------------

/// An `EdwardsPoint` represents a point on the Edwards form of Curve25519.
#[derive(Copy, Clone)]
#[allow(missing_docs)]
pub(crate) struct EdwardsPoint {
    pub(crate) X: FieldElement,
    pub(crate) Y: FieldElement,
    pub(crate) Z: FieldElement,
    pub(crate) T: FieldElement,
}

// ------------------------------------------------------------------------
// Constructors
// ------------------------------------------------------------------------

impl Identity for CompressedEdwardsY {
    fn identity() -> CompressedEdwardsY {
        CompressedEdwardsY([
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0,
        ])
    }
}

impl CompressedEdwardsY {
    /// Construct a `CompressedEdwardsY` from a slice of bytes.
    ///
    /// # Errors
    ///
    /// Returns [`TryFromSliceError`] if the input `bytes` slice does not have
    /// a length of 32.
    pub(crate) fn from_slice(bytes: &[u8]) -> Result<CompressedEdwardsY, TryFromSliceError> {
        bytes.try_into().map(CompressedEdwardsY)
    }
}

impl Identity for EdwardsPoint {
    fn identity() -> EdwardsPoint {
        EdwardsPoint {
            X: FieldElement::ZERO,
            Y: FieldElement::ONE,
            Z: FieldElement::ONE,
            T: FieldElement::ZERO,
        }
    }
}

// ------------------------------------------------------------------------
// Equality
// ------------------------------------------------------------------------

impl ConstantTimeEq for EdwardsPoint {
    fn ct_eq(&self, other: &EdwardsPoint) -> Choice {
        // We would like to check that the point (X/Z, Y/Z) is equal to
        // the point (X'/Z', Y'/Z') without converting into affine
        // coordinates (x, y) and (x', y'), which requires two inversions.
        // We have that X = xZ and X' = x'Z'. Thus, x = x' is equivalent to
        // (xZ)Z' = (x'Z')Z, and similarly for the y-coordinate.

        (&self.X * &other.Z).ct_eq(&(&other.X * &self.Z))
            & (&self.Y * &other.Z).ct_eq(&(&other.Y * &self.Z))
    }
}

impl PartialEq for EdwardsPoint {
    fn eq(&self, other: &EdwardsPoint) -> bool {
        self.ct_eq(other).into()
    }
}

impl Eq for EdwardsPoint {}

// ------------------------------------------------------------------------
// Point conversions
// ------------------------------------------------------------------------

impl EdwardsPoint {
    /// Convert to a ProjectiveNielsPoint
    pub(crate) fn as_projective_niels(&self) -> ProjectiveNielsPoint {
        ProjectiveNielsPoint {
            Y_plus_X: &self.Y + &self.X,
            Y_minus_X: &self.Y - &self.X,
            Z: self.Z,
            T2d: &self.T * &constants::EDWARDS_D2,
        }
    }

    /// Convert the representation of this point from extended
    /// coordinates to projective coordinates.
    ///
    /// Free.
    pub(crate) const fn as_projective(&self) -> ProjectivePoint {
        ProjectivePoint {
            X: self.X,
            Y: self.Y,
            Z: self.Z,
        }
    }

    /// Dehomogenize to a AffineNielsPoint.
    /// Mainly for testing.
    pub(crate) fn as_affine_niels(&self) -> AffineNielsPoint {
        let recip = self.Z.invert();
        let x = &self.X * &recip;
        let y = &self.Y * &recip;
        let xy2d = &(&x * &y) * &constants::EDWARDS_D2;
        AffineNielsPoint {
            y_plus_x: &y + &x,
            y_minus_x: &y - &x,
            xy2d,
        }
    }

    #[cfg(feature = "digest")]
    /// Maps the digest of the input bytes to the curve. This is NOT a hash-to-curve function, as
    /// it produces points with a non-uniform distribution. Rather, it performs something that
    /// resembles (but is not) half of the
    /// [`hash_to_curve`](https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-16.html#section-3-4.2.1)
    /// function from the Elligator2 spec.
    #[deprecated(
        since = "4.0.0",
        note = "previously named `hash_from_bytes`, this is not a secure hash function"
    )]
    pub(crate) fn nonspec_map_to_curve<D>(bytes: &[u8]) -> EdwardsPoint
    where
        D: Digest<OutputSize = U64> + Default,
    {
        let mut hash = D::new();
        hash.update(bytes);
        let h = hash.finalize();
        let mut res = [0u8; 32];
        res.copy_from_slice(&h[..32]);

        let sign_bit = (res[31] & 0x80) >> 7;

        let fe = FieldElement::from_bytes(&res);

        let M1 = crate::montgomery::elligator_encode(&fe);
        let E1_opt = M1.to_edwards(sign_bit);

        E1_opt
            .expect("Montgomery conversion to Edwards point in Elligator failed")
            .mul_by_cofactor()
    }
}

// ------------------------------------------------------------------------
// Doubling
// ------------------------------------------------------------------------

impl EdwardsPoint {
    /// Add this point to itself.
    pub(crate) fn double(&self) -> EdwardsPoint {
        self.as_projective().double().as_extended()
    }
}

// ------------------------------------------------------------------------
// Addition and Subtraction
// ------------------------------------------------------------------------

impl<'a, 'b> Add<&'b EdwardsPoint> for &'a EdwardsPoint {
    type Output = EdwardsPoint;
    fn add(self, other: &'b EdwardsPoint) -> EdwardsPoint {
        (self + &other.as_projective_niels()).as_extended()
    }
}

define_add_variants!(
    LHS = EdwardsPoint,
    RHS = EdwardsPoint,
    Output = EdwardsPoint
);

impl<'b> AddAssign<&'b EdwardsPoint> for EdwardsPoint {
    fn add_assign(&mut self, _rhs: &'b EdwardsPoint) {
        *self = (self as &EdwardsPoint) + _rhs;
    }
}

define_add_assign_variants!(LHS = EdwardsPoint, RHS = EdwardsPoint);

impl<'a, 'b> Sub<&'b EdwardsPoint> for &'a EdwardsPoint {
    type Output = EdwardsPoint;
    fn sub(self, other: &'b EdwardsPoint) -> EdwardsPoint {
        (self - &other.as_projective_niels()).as_extended()
    }
}

define_sub_variants!(
    LHS = EdwardsPoint,
    RHS = EdwardsPoint,
    Output = EdwardsPoint
);

impl<'b> SubAssign<&'b EdwardsPoint> for EdwardsPoint {
    fn sub_assign(&mut self, _rhs: &'b EdwardsPoint) {
        *self = (self as &EdwardsPoint) - _rhs;
    }
}

define_sub_assign_variants!(LHS = EdwardsPoint, RHS = EdwardsPoint);

impl<T> Sum<T> for EdwardsPoint
where
    T: Borrow<EdwardsPoint>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(EdwardsPoint::identity(), |acc, item| acc + item.borrow())
    }
}

// ------------------------------------------------------------------------
// Negation
// ------------------------------------------------------------------------

impl<'a> Neg for &'a EdwardsPoint {
    type Output = EdwardsPoint;

    fn neg(self) -> EdwardsPoint {
        EdwardsPoint {
            X: -(&self.X),
            Y: self.Y,
            Z: self.Z,
            T: -(&self.T),
        }
    }
}

impl Neg for EdwardsPoint {
    type Output = EdwardsPoint;

    fn neg(self) -> EdwardsPoint {
        -&self
    }
}

// ------------------------------------------------------------------------
// Scalar multiplication
// ------------------------------------------------------------------------

impl<'b> MulAssign<&'b Scalar> for EdwardsPoint {
    fn mul_assign(&mut self, scalar: &'b Scalar) {
        let result = (self as &EdwardsPoint) * scalar;
        *self = result;
    }
}

define_mul_assign_variants!(LHS = EdwardsPoint, RHS = Scalar);

define_mul_variants!(LHS = EdwardsPoint, RHS = Scalar, Output = EdwardsPoint);
define_mul_variants!(LHS = Scalar, RHS = EdwardsPoint, Output = EdwardsPoint);

impl<'a, 'b> Mul<&'b Scalar> for &'a EdwardsPoint {
    type Output = EdwardsPoint;
    /// Scalar multiplication: compute `scalar * self`.
    ///
    /// For scalar multiplication of a basepoint,
    /// `EdwardsBasepointTable` is approximately 4x faster.
    fn mul(self, scalar: &'b Scalar) -> EdwardsPoint {
        crate::backend::variable_base_mul(self, scalar)
    }
}

impl<'a, 'b> Mul<&'b EdwardsPoint> for &'a Scalar {
    type Output = EdwardsPoint;

    /// Scalar multiplication: compute `scalar * self`.
    ///
    /// For scalar multiplication of a basepoint,
    /// `EdwardsBasepointTable` is approximately 4x faster.
    fn mul(self, point: &'b EdwardsPoint) -> EdwardsPoint {
        point * self
    }
}





#[cfg(feature = "precomputed-tables")]
macro_rules! impl_basepoint_table {
    (Name = $name:ident, LookupTable = $table:ident, Point = $point:ty, Radix = $radix:expr, Additions = $adds:expr) => {
        /// A precomputed table of multiples of a basepoint, for accelerating
        /// fixed-base scalar multiplication.  One table, for the Ed25519
        /// basepoint, is provided in the [`constants`] module.
        ///
        /// The basepoint tables are reasonably large, so they should probably be boxed.
        ///
        /// The sizes for the tables and the number of additions required for one scalar
        /// multiplication are as follows:
        ///
        /// * [`EdwardsBasepointTableRadix16`]: 30KB, 64A
        ///   (this is the default size, and is used for
        ///   [`constants::ED25519_BASEPOINT_TABLE`])
        /// * [`EdwardsBasepointTableRadix64`]: 120KB, 43A
        /// * [`EdwardsBasepointTableRadix128`]: 240KB, 37A
        /// * [`EdwardsBasepointTableRadix256`]: 480KB, 33A
        ///
        /// # Why 33 additions for radix-256?
        ///
        /// Normally, the radix-256 tables would allow for only 32 additions per scalar
        /// multiplication.  However, due to the fact that standardised definitions of
        /// legacy protocols—such as x25519—require allowing unreduced 255-bit scalars
        /// invariants, when converting such an unreduced scalar's representation to
        /// radix-\\(2^{8}\\), we cannot guarantee the carry bit will fit in the last
        /// coefficient (the coefficients are `i8`s).  When, \\(w\\), the power-of-2 of
        /// the radix, is \\(w < 8\\), we can fold the final carry onto the last
        /// coefficient, \\(d\\), because \\(d < 2^{w/2}\\), so
        /// $$
        ///     d + carry \cdot 2^{w} = d + 1 \cdot 2^{w} < 2^{w+1} < 2^{8}
        /// $$
        /// When \\(w = 8\\), we can't fit \\(carry \cdot 2^{w}\\) into an `i8`, so we
        /// add the carry bit onto an additional coefficient.
        #[derive(Clone)]
        #[repr(transparent)]
        pub(crate) struct $name(pub(crate) [$table<AffineNielsPoint>; 32]);

        impl BasepointTable for $name {
            type Point = $point;

            /// The computation uses Pippeneger's algorithm, as described for the
            /// specific case of radix-16 on page 13 of the Ed25519 paper.
            ///
            /// # Piggenger's Algorithm Generalised
            ///
            /// Write the scalar \\(a\\) in radix-\\(w\\), where \\(w\\) is a power of
            /// 2, with coefficients in \\([\frac{-w}{2},\frac{w}{2})\\), i.e.,
            /// $$
            ///     a = a\_0 + a\_1 w\^1 + \cdots + a\_{x} w\^{x},
            /// $$
            /// with
            /// $$
            /// \begin{aligned}
            ///     \frac{-w}{2} \leq a_i < \frac{w}{2}
            ///     &&\cdots&&
            ///     \frac{-w}{2} \leq a\_{x} \leq \frac{w}{2}
            /// \end{aligned}
            /// $$
            /// and the number of additions, \\(x\\), is given by
            /// \\(x = \lceil \frac{256}{w} \rceil\\). Then
            /// $$
            ///     a B = a\_0 B + a\_1 w\^1 B + \cdots + a\_{x-1} w\^{x-1} B.
            /// $$
            /// Grouping even and odd coefficients gives
            /// $$
            /// \begin{aligned}
            ///     a B = \quad a\_0 w\^0 B +& a\_2 w\^2 B + \cdots + a\_{x-2} w\^{x-2} B    \\\\
            ///               + a\_1 w\^1 B +& a\_3 w\^3 B + \cdots + a\_{x-1} w\^{x-1} B    \\\\
            ///         = \quad(a\_0 w\^0 B +& a\_2 w\^2 B + \cdots + a\_{x-2} w\^{x-2} B)   \\\\
            ///             + w(a\_1 w\^0 B +& a\_3 w\^2 B + \cdots + a\_{x-1} w\^{x-2} B).  \\\\
            /// \end{aligned}
            /// $$
            /// For each \\(i = 0 \ldots 31\\), we create a lookup table of
            /// $$
            /// [w\^{2i} B, \ldots, \frac{w}{2}\cdot w\^{2i} B],
            /// $$
            /// and use it to select \\( y \cdot w\^{2i} \cdot B \\) in constant time.
            ///
            /// The radix-\\(w\\) representation requires that the scalar is bounded
            /// by \\(2\^{255}\\), which is always the case.
            ///
            /// The above algorithm is trivially generalised to other powers-of-2 radices.
            fn mul_base(&self, scalar: &Scalar) -> $point {
                let a = scalar.as_radix_2w($radix);

                let tables = &self.0;
                let mut P = <$point>::identity();

                for i in (0..$adds).filter(|x| x % 2 == 1) {
                    P = (&P + &tables[i / 2].select(a[i])).as_extended();
                }

                P = P.mul_by_pow_2($radix);

                for i in (0..$adds).filter(|x| x % 2 == 0) {
                    P = (&P + &tables[i / 2].select(a[i])).as_extended();
                }

                P
            }
        }

        impl<'a, 'b> Mul<&'b Scalar> for &'a $name {
            type Output = $point;

            /// Construct an `EdwardsPoint` from a `Scalar` \\(a\\) by
            /// computing the multiple \\(aB\\) of this basepoint \\(B\\).
            fn mul(self, scalar: &'b Scalar) -> $point {
                // delegate to a private function so that its documentation appears in internal docs
                self.mul_base(scalar)
            }
        }

        impl<'a, 'b> Mul<&'a $name> for &'b Scalar {
            type Output = $point;

            /// Construct an `EdwardsPoint` from a `Scalar` \\(a\\) by
            /// computing the multiple \\(aB\\) of this basepoint \\(B\\).
            fn mul(self, basepoint_table: &'a $name) -> $point {
                basepoint_table * self
            }
        }

        impl Debug for $name {
            fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
                write!(f, "{:?}([\n", stringify!($name))?;
                for i in 0..32 {
                    write!(f, "\t{:?},\n", &self.0[i])?;
                }
                write!(f, "])")
            }
        }
    };
} // End macro_rules! impl_basepoint_table

// The number of additions required is ceil(256/w) where w is the radix representation.
cfg_if! {
    if #[cfg(feature = "precomputed-tables")] {
        impl_basepoint_table! {
            Name = EdwardsBasepointTable,
            LookupTable = LookupTableRadix16,
            Point = EdwardsPoint,
            Radix = 4,
            Additions = 64
        }
        impl_basepoint_table! {
            Name = EdwardsBasepointTableRadix32,
            LookupTable = LookupTableRadix32,
            Point = EdwardsPoint,
            Radix = 5,
            Additions = 52
        }
        impl_basepoint_table! {
            Name = EdwardsBasepointTableRadix64,
            LookupTable = LookupTableRadix64,
            Point = EdwardsPoint,
            Radix = 6,
            Additions = 43
        }
        impl_basepoint_table! {
            Name = EdwardsBasepointTableRadix128,
            LookupTable = LookupTableRadix128,
            Point = EdwardsPoint,
            Radix = 7,
            Additions = 37
        }
        impl_basepoint_table! {
            Name = EdwardsBasepointTableRadix256,
            LookupTable = LookupTableRadix256,
            Point = EdwardsPoint,
            Radix = 8,
            Additions = 33
        }

    }
}


impl EdwardsPoint {
    /// Compute \\([2\^k] P \\) by successive doublings. Requires \\( k > 0 \\).
    pub(crate) fn mul_by_pow_2(&self, k: u32) -> EdwardsPoint {
        debug_assert!(k > 0);
        let mut r: CompletedPoint;
        let mut s = self.as_projective();
        for _ in 0..(k - 1) {
            r = s.double();
            s = r.as_projective();
        }
        // Unroll last iteration so we can go directly as_extended()
        s.double().as_extended()
    }
}

// ------------------------------------------------------------------------
// Debug traits
// ------------------------------------------------------------------------

impl Debug for EdwardsPoint {
    fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
        write!(
            f,
            "EdwardsPoint{{\n\tX: {:?},\n\tY: {:?},\n\tZ: {:?},\n\tT: {:?}\n}}",
            &self.X, &self.Y, &self.Z, &self.T
        )
    }
}
