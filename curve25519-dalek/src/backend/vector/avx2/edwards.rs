// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - isis agora lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>

//! Parallel Edwards Arithmetic for Curve25519.
//!
//! This module currently has two point types:
//!
//! * `ExtendedPoint`: a point stored in vector-friendly format, with
//! vectorized doubling and addition;
//!
//! * `CachedPoint`: used for readdition.
//!
//! Details on the formulas can be found in the documentation for the
//! parent `avx2` module.
//!
//! This API is designed to be safe: vectorized points can only be
//! created from serial points (which do validation on decompression),
//! and operations on valid points return valid points, so invalid
//! point states should be unrepresentable.
//!
//! This design goal is met, with one exception: the `Neg`
//! implementation for the `CachedPoint` performs a lazy negation, so
//! that subtraction can be efficiently implemented as a negation and
//! an addition.  Repeatedly negating a `CachedPoint` will cause its
//! coefficients to grow and eventually overflow.  Repeatedly negating
//! a point should not be necessary anyways.

#![allow(non_snake_case)]


use curve25519_dalek_derive::unsafe_target_feature;


use super::field::{FieldElement2625x4, Lanes, Shuffle};

/// A point on Curve25519, using parallel Edwards formulas for curve
/// operations.
///
/// # Invariant
///
/// The coefficients of an `ExtendedPoint` are bounded with
/// \\( b < 0.007 \\).
#[derive(Copy, Clone, Debug)]
pub(crate) struct ExtendedPoint(pub(super) FieldElement2625x4);



#[unsafe_target_feature("avx2")]
impl ExtendedPoint {
    /// Compute the double of this point.
    pub(crate) fn double(&self) -> ExtendedPoint {
        // Want to compute (X1 Y1 Z1 X1+Y1).
        // Not sure how to do this less expensively than computing
        // (X1 Y1 Z1 T1) --(256bit shuffle)--> (X1 Y1 X1 Y1)
        // (X1 Y1 X1 Y1) --(2x128b shuffle)--> (Y1 X1 Y1 X1)
        // and then adding.

        // Set tmp0 = (X1 Y1 X1 Y1)
        let mut tmp0 = self.0.shuffle(Shuffle::ABAB);

        // Set tmp1 = (Y1 X1 Y1 X1)
        let mut tmp1 = tmp0.shuffle(Shuffle::BADC);

        // Set tmp0 = (X1 Y1 Z1 X1+Y1)
        tmp0 = self.0.blend(tmp0 + tmp1, Lanes::D);

        // Set tmp1 = tmp0^2, negating the D values
        tmp1 = tmp0.square_and_negate_D();
        // Now tmp1 = (S1 S2 S3 -S4) with b < 0.007

        // See discussion of bounds in the module-level documentation.
        // We want to compute
        //
        //    + | S1 | S1 | S1 | S1 |
        //    + | S2 |    |    | S2 |
        //    + |    |    | S3 |    |
        //    + |    |    | S3 |    |
        //    + |    |    |    |-S4 |
        //    + |    | 2p | 2p |    |
        //    - |    | S2 | S2 |    |
        //    =======================
        //        S5   S6   S8   S9

        let zero = FieldElement2625x4::ZERO;
        let S_1 = tmp1.shuffle(Shuffle::AAAA);
        let S_2 = tmp1.shuffle(Shuffle::BBBB);

        tmp0 = zero.blend(tmp1 + tmp1, Lanes::C);
        // tmp0 = (0, 0,  2S_3, 0)
        tmp0 = tmp0.blend(tmp1, Lanes::D);
        // tmp0 = (0, 0,  2S_3, -S_4)
        tmp0 = tmp0 + S_1;
        // tmp0 = (  S_1,   S_1, S_1 + 2S_3, S_1 - S_4)
        tmp0 = tmp0 + zero.blend(S_2, Lanes::AD);
        // tmp0 = (S_1 + S_2,   S_1, S_1 + 2S_3, S_1 + S_2 - S_4)
        tmp0 = tmp0 + zero.blend(S_2.negate_lazy(), Lanes::BC);
        // tmp0 = (S_1 + S_2, S_1 - S_2, S_1 - S_2 + 2S_3, S_1 + S_2 - S_4)
        //    b < (     1.01,       1.6,             2.33,             1.6)
        // Now tmp0 = (S_5, S_6, S_8, S_9)

        // Set tmp1 = ( S_9,  S_6,  S_6,  S_9)
        //        b < ( 1.6,  1.6,  1.6,  1.6)
        tmp1 = tmp0.shuffle(Shuffle::DBBD);
        // Set tmp0 = ( S_8,  S_5,  S_8,  S_5)
        //        b < (2.33, 1.01, 2.33, 1.01)
        tmp0 = tmp0.shuffle(Shuffle::CACA);

        // Bounds on (tmp0, tmp1) are (2.33, 1.6) < (2.5, 1.75).
        ExtendedPoint(&tmp0 * &tmp1)
    }

    pub(crate) fn mul_by_pow_2(&self, k: u32) -> ExtendedPoint {
        let mut tmp: ExtendedPoint = *self;
        for _ in 0..k {
            tmp = tmp.double();
        }
        tmp
    }
}

/// A cached point with some precomputed variables used for readdition.
///
/// # Warning
///
/// It is not safe to negate this point more than once.
///
/// # Invariant
///
/// As long as the `CachedPoint` is not repeatedly negated, its
/// coefficients will be bounded with \\( b < 1.0 \\).
#[derive(Copy, Clone, Debug)]
pub(crate) struct CachedPoint(pub(super) FieldElement2625x4);
