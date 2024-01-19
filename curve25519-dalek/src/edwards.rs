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

use subtle::Choice;

use crate::constants;

use crate::field::FieldElement;


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





impl CompressedEdwardsY {
    /// View this `CompressedEdwardsY` as an array of bytes.
    pub(crate) const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }

    pub fn is_curve_point(&self) -> bool {
        let is_valid_y_coord = decompress::step_1(self).0;
        is_valid_y_coord.into()
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
}


// ------------------------------------------------------------------------
// Constructors
// ------------------------------------------------------------------------



impl CompressedEdwardsY {
    /// Construct a `CompressedEdwardsY` from a slice of bytes.
    ///
    /// # Errors
    ///
    /// Returns [`TryFromSliceError`] if the input `bytes` slice does not have
    /// a length of 32.
    pub fn from_slice(bytes: &[u8]) -> Result<CompressedEdwardsY, TryFromSliceError> {
        bytes.try_into().map(CompressedEdwardsY)
    }
}
