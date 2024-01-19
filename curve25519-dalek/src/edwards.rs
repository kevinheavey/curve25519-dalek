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
use subtle::ConditionallySelectable;
use subtle::ConstantTimeEq;

#[cfg(feature = "zeroize")]
use zeroize::Zeroize;

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

use crate::traits::ValidityCheck;
use crate::traits::{Identity};

#[cfg(feature = "alloc")]
use crate::traits::MultiscalarMul;
#[cfg(feature = "alloc")]
use crate::traits::{VartimeMultiscalarMul, VartimePrecomputedMultiscalarMul};

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
// Serde support
// ------------------------------------------------------------------------
// Serializes to and from `EdwardsPoint` directly, doing compression
// and decompression internally.  This means that users can create
// structs containing `EdwardsPoint`s and use Serde's derived
// serializers to serialize those structures.

#[cfg(feature = "serde")]
use serde::de::Visitor;
#[cfg(feature = "serde")]
use serde::{self, Deserialize, Deserializer, Serialize, Serializer};

#[cfg(feature = "serde")]
impl Serialize for EdwardsPoint {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(32)?;
        for byte in self.compress().as_bytes().iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serde")]
impl Serialize for CompressedEdwardsY {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(32)?;
        for byte in self.as_bytes().iter() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for EdwardsPoint {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct EdwardsPointVisitor;

        impl<'de> Visitor<'de> for EdwardsPointVisitor {
            type Value = EdwardsPoint;

            fn expecting(&self, formatter: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
                formatter.write_str("a valid point in Edwards y + sign format")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<EdwardsPoint, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 32];
                #[allow(clippy::needless_range_loop)]
                for i in 0..32 {
                    bytes[i] = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 32 bytes"))?;
                }
                CompressedEdwardsY(bytes)
                    .decompress()
                    .ok_or_else(|| serde::de::Error::custom("decompression failed"))
            }
        }

        deserializer.deserialize_tuple(32, EdwardsPointVisitor)
    }
}

#[cfg(feature = "serde")]
impl<'de> Deserialize<'de> for CompressedEdwardsY {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct CompressedEdwardsYVisitor;

        impl<'de> Visitor<'de> for CompressedEdwardsYVisitor {
            type Value = CompressedEdwardsY;

            fn expecting(&self, formatter: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
                formatter.write_str("32 bytes of data")
            }

            fn visit_seq<A>(self, mut seq: A) -> Result<CompressedEdwardsY, A::Error>
            where
                A: serde::de::SeqAccess<'de>,
            {
                let mut bytes = [0u8; 32];
                #[allow(clippy::needless_range_loop)]
                for i in 0..32 {
                    bytes[i] = seq
                        .next_element()?
                        .ok_or_else(|| serde::de::Error::invalid_length(i, &"expected 32 bytes"))?;
                }
                Ok(CompressedEdwardsY(bytes))
            }
        }

        deserializer.deserialize_tuple(32, CompressedEdwardsYVisitor)
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

impl Default for CompressedEdwardsY {
    fn default() -> CompressedEdwardsY {
        CompressedEdwardsY::identity()
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

impl Default for EdwardsPoint {
    fn default() -> EdwardsPoint {
        EdwardsPoint::identity()
    }
}

// ------------------------------------------------------------------------
// Zeroize implementations for wiping points from memory
// ------------------------------------------------------------------------

#[cfg(feature = "zeroize")]
impl Zeroize for CompressedEdwardsY {
    /// Reset this `CompressedEdwardsY` to the compressed form of the identity element.
    fn zeroize(&mut self) {
        self.0.zeroize();
        self.0[0] = 1;
    }
}

#[cfg(feature = "zeroize")]
impl Zeroize for EdwardsPoint {
    /// Reset this `CompressedEdwardsPoint` to the identity element.
    fn zeroize(&mut self) {
        self.X.zeroize();
        self.Y = FieldElement::ONE;
        self.Z = FieldElement::ONE;
        self.T.zeroize();
    }
}

// ------------------------------------------------------------------------
// Validity checks (for debugging, not CT)
// ------------------------------------------------------------------------

impl ValidityCheck for EdwardsPoint {
    fn is_valid(&self) -> bool {
        let point_on_curve = self.as_projective().is_valid();
        let on_segre_image = (&self.X * &self.Y) == (&self.Z * &self.T);

        point_on_curve && on_segre_image
    }
}

// ------------------------------------------------------------------------
// Constant-time assignment
// ------------------------------------------------------------------------

impl ConditionallySelectable for EdwardsPoint {
    fn conditional_select(a: &EdwardsPoint, b: &EdwardsPoint, choice: Choice) -> EdwardsPoint {
        EdwardsPoint {
            X: FieldElement::conditional_select(&a.X, &b.X, choice),
            Y: FieldElement::conditional_select(&a.Y, &b.Y, choice),
            Z: FieldElement::conditional_select(&a.Z, &b.Z, choice),
            T: FieldElement::conditional_select(&a.T, &b.T, choice),
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


// ------------------------------------------------------------------------
// Multiscalar Multiplication impls
// ------------------------------------------------------------------------

// These use the iterator's size hint and the target settings to
// forward to a specific backend implementation.

#[cfg(feature = "alloc")]
impl MultiscalarMul for EdwardsPoint {
    type Point = EdwardsPoint;

    fn multiscalar_mul<I, J>(scalars: I, points: J) -> EdwardsPoint
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
        J: IntoIterator,
        J::Item: Borrow<EdwardsPoint>,
    {
        // Sanity-check lengths of input iterators
        let mut scalars = scalars.into_iter();
        let mut points = points.into_iter();

        // Lower and upper bounds on iterators
        let (s_lo, s_hi) = scalars.by_ref().size_hint();
        let (p_lo, p_hi) = points.by_ref().size_hint();

        // They should all be equal
        assert_eq!(s_lo, p_lo);
        assert_eq!(s_hi, Some(s_lo));
        assert_eq!(p_hi, Some(p_lo));

        // Now we know there's a single size.  When we do
        // size-dependent algorithm dispatch, use this as the hint.
        let _size = s_lo;

        crate::backend::straus_multiscalar_mul(scalars, points)
    }
}

#[cfg(feature = "alloc")]
impl VartimeMultiscalarMul for EdwardsPoint {
    type Point = EdwardsPoint;

    fn optional_multiscalar_mul<I, J>(scalars: I, points: J) -> Option<EdwardsPoint>
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
        J: IntoIterator<Item = Option<EdwardsPoint>>,
    {
        // Sanity-check lengths of input iterators
        let mut scalars = scalars.into_iter();
        let mut points = points.into_iter();

        // Lower and upper bounds on iterators
        let (s_lo, s_hi) = scalars.by_ref().size_hint();
        let (p_lo, p_hi) = points.by_ref().size_hint();

        // They should all be equal
        assert_eq!(s_lo, p_lo);
        assert_eq!(s_hi, Some(s_lo));
        assert_eq!(p_hi, Some(p_lo));

        // Now we know there's a single size.
        // Use this as the hint to decide which algorithm to use.
        let size = s_lo;

        if size < 190 {
            crate::backend::straus_optional_multiscalar_mul(scalars, points)
        } else {
            crate::backend::pippenger_optional_multiscalar_mul(scalars, points)
        }
    }
}

/// Precomputation for variable-time multiscalar multiplication with `EdwardsPoint`s.
// This wraps the inner implementation in a facade type so that we can
// decouple stability of the inner type from the stability of the
// outer type.
#[cfg(feature = "alloc")]
pub(crate) struct VartimeEdwardsPrecomputation(crate::backend::VartimePrecomputedStraus);

#[cfg(feature = "alloc")]
impl VartimePrecomputedMultiscalarMul for VartimeEdwardsPrecomputation {
    type Point = EdwardsPoint;

    fn new<I>(static_points: I) -> Self
    where
        I: IntoIterator,
        I::Item: Borrow<Self::Point>,
    {
        Self(crate::backend::VartimePrecomputedStraus::new(static_points))
    }

    fn optional_mixed_multiscalar_mul<I, J, K>(
        &self,
        static_scalars: I,
        dynamic_scalars: J,
        dynamic_points: K,
    ) -> Option<Self::Point>
    where
        I: IntoIterator,
        I::Item: Borrow<Scalar>,
        J: IntoIterator,
        J::Item: Borrow<Scalar>,
        K: IntoIterator<Item = Option<Self::Point>>,
    {
        self.0
            .optional_mixed_multiscalar_mul(static_scalars, dynamic_scalars, dynamic_points)
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

            /// Create a table of precomputed multiples of `basepoint`.
            fn create(basepoint: &$point) -> $name {
                // XXX use init_with
                let mut table = $name([$table::default(); 32]);
                let mut P = *basepoint;
                for i in 0..32 {
                    // P = (2w)^i * B
                    table.0[i] = $table::from(&P);
                    P = P.mul_by_pow_2($radix + $radix);
                }
                table
            }

            /// Get the basepoint for this table as an `EdwardsPoint`.
            fn basepoint(&self) -> $point {
                // self.0[0].select(1) = 1*(16^2)^0*B
                // but as an `AffineNielsPoint`, so add identity to convert to extended.
                (&<$point>::identity() + &self.0[0].select(1)).as_extended()
            }

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

        /// A type-alias for [`EdwardsBasepointTable`] because the latter is
        /// used as a constructor in the [`constants`] module.
        //
        // Same as for `LookupTableRadix16`, we have to define `EdwardsBasepointTable`
        // first, because it's used as a constructor, and then provide a type alias for
        // it.
        pub(crate) type EdwardsBasepointTableRadix16 = EdwardsBasepointTable;
    }
}

#[cfg(feature = "precomputed-tables")]
macro_rules! impl_basepoint_table_conversions {
    (LHS = $lhs:ty, RHS = $rhs:ty) => {
        impl<'a> From<&'a $lhs> for $rhs {
            fn from(table: &'a $lhs) -> $rhs {
                <$rhs>::create(&table.basepoint())
            }
        }

        impl<'a> From<&'a $rhs> for $lhs {
            fn from(table: &'a $rhs) -> $lhs {
                <$lhs>::create(&table.basepoint())
            }
        }
    };
}

cfg_if! {
    if #[cfg(feature = "precomputed-tables")] {
        // Conversions from radix 16
        impl_basepoint_table_conversions! {
            LHS = EdwardsBasepointTableRadix16,
            RHS = EdwardsBasepointTableRadix32
        }
        impl_basepoint_table_conversions! {
            LHS = EdwardsBasepointTableRadix16,
            RHS = EdwardsBasepointTableRadix64
        }
        impl_basepoint_table_conversions! {
            LHS = EdwardsBasepointTableRadix16,
            RHS = EdwardsBasepointTableRadix128
        }
        impl_basepoint_table_conversions! {
            LHS = EdwardsBasepointTableRadix16,
            RHS = EdwardsBasepointTableRadix256
        }

        // Conversions from radix 32
        impl_basepoint_table_conversions! {
            LHS = EdwardsBasepointTableRadix32,
            RHS = EdwardsBasepointTableRadix64
        }
        impl_basepoint_table_conversions! {
            LHS = EdwardsBasepointTableRadix32,
            RHS = EdwardsBasepointTableRadix128
        }
        impl_basepoint_table_conversions! {
            LHS = EdwardsBasepointTableRadix32,
            RHS = EdwardsBasepointTableRadix256
        }

        // Conversions from radix 64
        impl_basepoint_table_conversions! {
            LHS = EdwardsBasepointTableRadix64,
            RHS = EdwardsBasepointTableRadix128
        }
        impl_basepoint_table_conversions! {
            LHS = EdwardsBasepointTableRadix64,
            RHS = EdwardsBasepointTableRadix256
        }

        // Conversions from radix 128
        impl_basepoint_table_conversions! {
            LHS = EdwardsBasepointTableRadix128,
            RHS = EdwardsBasepointTableRadix256
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

// ------------------------------------------------------------------------
// group traits
// ------------------------------------------------------------------------

// Use the full trait path to avoid Group::identity overlapping Identity::identity in the
// rest of the module (e.g. tests).
#[cfg(feature = "group")]
impl group::Group for EdwardsPoint {
    type Scalar = Scalar;

    fn random(mut rng: impl RngCore) -> Self {
        let mut repr = CompressedEdwardsY([0u8; 32]);
        loop {
            rng.fill_bytes(&mut repr.0);
            if let Some(p) = repr.decompress() {
                if !IsIdentity::is_identity(&p) {
                    break p;
                }
            }
        }
    }

    fn identity() -> Self {
        Identity::identity()
    }

    fn generator() -> Self {
        constants::ED25519_BASEPOINT_POINT
    }

    fn is_identity(&self) -> Choice {
        self.ct_eq(&Identity::identity())
    }

    fn double(&self) -> Self {
        self.double()
    }
}

#[cfg(feature = "group")]
impl GroupEncoding for EdwardsPoint {
    type Repr = [u8; 32];

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        let repr = CompressedEdwardsY(*bytes);
        let (is_valid_y_coord, X, Y, Z) = decompress::step_1(&repr);
        CtOption::new(decompress::step_2(&repr, X, Y, Z), is_valid_y_coord)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        // Just use the checked API; there are no checks we can skip.
        Self::from_bytes(bytes)
    }

    fn to_bytes(&self) -> Self::Repr {
        self.compress().to_bytes()
    }
}

/// A `SubgroupPoint` represents a point on the Edwards form of Curve25519, that is
/// guaranteed to be in the prime-order subgroup.
#[cfg(feature = "group")]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SubgroupPoint(EdwardsPoint);

#[cfg(feature = "group")]
impl From<SubgroupPoint> for EdwardsPoint {
    fn from(p: SubgroupPoint) -> Self {
        p.0
    }
}

#[cfg(feature = "group")]
impl Neg for SubgroupPoint {
    type Output = Self;

    fn neg(self) -> Self::Output {
        SubgroupPoint(-self.0)
    }
}

#[cfg(feature = "group")]
impl Add<&SubgroupPoint> for &SubgroupPoint {
    type Output = SubgroupPoint;
    fn add(self, other: &SubgroupPoint) -> SubgroupPoint {
        SubgroupPoint(self.0 + other.0)
    }
}

#[cfg(feature = "group")]
define_add_variants!(
    LHS = SubgroupPoint,
    RHS = SubgroupPoint,
    Output = SubgroupPoint
);

#[cfg(feature = "group")]
impl Add<&SubgroupPoint> for &EdwardsPoint {
    type Output = EdwardsPoint;
    fn add(self, other: &SubgroupPoint) -> EdwardsPoint {
        self + other.0
    }
}

#[cfg(feature = "group")]
define_add_variants!(
    LHS = EdwardsPoint,
    RHS = SubgroupPoint,
    Output = EdwardsPoint
);

#[cfg(feature = "group")]
impl AddAssign<&SubgroupPoint> for SubgroupPoint {
    fn add_assign(&mut self, rhs: &SubgroupPoint) {
        self.0 += rhs.0
    }
}

#[cfg(feature = "group")]
define_add_assign_variants!(LHS = SubgroupPoint, RHS = SubgroupPoint);

#[cfg(feature = "group")]
impl AddAssign<&SubgroupPoint> for EdwardsPoint {
    fn add_assign(&mut self, rhs: &SubgroupPoint) {
        *self += rhs.0
    }
}

#[cfg(feature = "group")]
define_add_assign_variants!(LHS = EdwardsPoint, RHS = SubgroupPoint);

#[cfg(feature = "group")]
impl Sub<&SubgroupPoint> for &SubgroupPoint {
    type Output = SubgroupPoint;
    fn sub(self, other: &SubgroupPoint) -> SubgroupPoint {
        SubgroupPoint(self.0 - other.0)
    }
}

#[cfg(feature = "group")]
define_sub_variants!(
    LHS = SubgroupPoint,
    RHS = SubgroupPoint,
    Output = SubgroupPoint
);

#[cfg(feature = "group")]
impl Sub<&SubgroupPoint> for &EdwardsPoint {
    type Output = EdwardsPoint;
    fn sub(self, other: &SubgroupPoint) -> EdwardsPoint {
        self - other.0
    }
}

#[cfg(feature = "group")]
define_sub_variants!(
    LHS = EdwardsPoint,
    RHS = SubgroupPoint,
    Output = EdwardsPoint
);

#[cfg(feature = "group")]
impl SubAssign<&SubgroupPoint> for SubgroupPoint {
    fn sub_assign(&mut self, rhs: &SubgroupPoint) {
        self.0 -= rhs.0;
    }
}

#[cfg(feature = "group")]
define_sub_assign_variants!(LHS = SubgroupPoint, RHS = SubgroupPoint);

#[cfg(feature = "group")]
impl SubAssign<&SubgroupPoint> for EdwardsPoint {
    fn sub_assign(&mut self, rhs: &SubgroupPoint) {
        *self -= rhs.0;
    }
}

#[cfg(feature = "group")]
define_sub_assign_variants!(LHS = EdwardsPoint, RHS = SubgroupPoint);

#[cfg(feature = "group")]
impl<T> Sum<T> for SubgroupPoint
where
    T: Borrow<SubgroupPoint>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        use group::Group;
        iter.fold(SubgroupPoint::identity(), |acc, item| acc + item.borrow())
    }
}

#[cfg(feature = "group")]
impl Mul<&Scalar> for &SubgroupPoint {
    type Output = SubgroupPoint;

    /// Scalar multiplication: compute `scalar * self`.
    ///
    /// For scalar multiplication of a basepoint,
    /// `EdwardsBasepointTable` is approximately 4x faster.
    fn mul(self, scalar: &Scalar) -> SubgroupPoint {
        SubgroupPoint(self.0 * scalar)
    }
}

#[cfg(feature = "group")]
define_mul_variants!(LHS = Scalar, RHS = SubgroupPoint, Output = SubgroupPoint);

#[cfg(feature = "group")]
impl Mul<&SubgroupPoint> for &Scalar {
    type Output = SubgroupPoint;

    /// Scalar multiplication: compute `scalar * self`.
    ///
    /// For scalar multiplication of a basepoint,
    /// `EdwardsBasepointTable` is approximately 4x faster.
    fn mul(self, point: &SubgroupPoint) -> SubgroupPoint {
        point * self
    }
}

#[cfg(feature = "group")]
define_mul_variants!(LHS = SubgroupPoint, RHS = Scalar, Output = SubgroupPoint);

#[cfg(feature = "group")]
impl MulAssign<&Scalar> for SubgroupPoint {
    fn mul_assign(&mut self, scalar: &Scalar) {
        self.0 *= scalar;
    }
}

#[cfg(feature = "group")]
define_mul_assign_variants!(LHS = SubgroupPoint, RHS = Scalar);

#[cfg(feature = "group")]
impl group::Group for SubgroupPoint {
    type Scalar = Scalar;

    fn random(mut rng: impl RngCore) -> Self {
        use group::ff::Field;

        // This will almost never loop, but `Group::random` is documented as returning a
        // non-identity element.
        let s = loop {
            let s: Scalar = Field::random(&mut rng);
            if !s.is_zero_vartime() {
                break s;
            }
        };

        // This gives an element of the prime-order subgroup.
        Self::generator() * s
    }

    fn identity() -> Self {
        SubgroupPoint(Identity::identity())
    }

    fn generator() -> Self {
        SubgroupPoint(EdwardsPoint::generator())
    }

    fn is_identity(&self) -> Choice {
        self.0.ct_eq(&Identity::identity())
    }

    fn double(&self) -> Self {
        SubgroupPoint(self.0.double())
    }
}

#[cfg(feature = "group")]
impl GroupEncoding for SubgroupPoint {
    type Repr = <EdwardsPoint as GroupEncoding>::Repr;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        EdwardsPoint::from_bytes(bytes).and_then(|p| p.into_subgroup())
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        EdwardsPoint::from_bytes_unchecked(bytes).and_then(|p| p.into_subgroup())
    }

    fn to_bytes(&self) -> Self::Repr {
        self.0.compress().to_bytes()
    }
}

#[cfg(feature = "group")]
impl PrimeGroup for SubgroupPoint {}

/// Ristretto has a cofactor of 1.
#[cfg(feature = "group")]
impl CofactorGroup for EdwardsPoint {
    type Subgroup = SubgroupPoint;

    fn clear_cofactor(&self) -> Self::Subgroup {
        SubgroupPoint(self.mul_by_cofactor())
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        CtOption::new(SubgroupPoint(self), CofactorGroup::is_torsion_free(&self))
    }

    fn is_torsion_free(&self) -> Choice {
        (self * constants::BASEPOINT_ORDER_PRIVATE).ct_eq(&Self::identity())
    }
}
