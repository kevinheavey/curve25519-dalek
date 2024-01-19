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

//! An implementation of 4-way vectorized 32bit field arithmetic using
//! AVX2.
//!
//! The `FieldElement2625x4` struct provides a vector of four field
//! elements, implemented using AVX2 operations.  Its API is designed
//! to abstract away the platform-dependent details, so that point
//! arithmetic can be implemented only in terms of a vector of field
//! elements.
//!
//! At this level, the API is optimized for speed and not safety.  The
//! `FieldElement2625x4` does not always perform reductions.  The pre-
//! and post-conditions on the bounds of the coefficients are
//! documented for each method, but it is the caller's responsibility
//! to ensure that there are no overflows.

#![allow(non_snake_case)]

#[allow(unused)]
const A_LANES64: u8 = 0b00_00_00_11;
#[allow(unused)]
const B_LANES64: u8 = 0b00_00_11_00;
#[allow(unused)]
const C_LANES64: u8 = 0b00_11_00_00;
#[allow(unused)]
const D_LANES64: u8 = 0b11_00_00_00;

use crate::backend::vector::packed_simd::u32x8;

use crate::backend::serial::u64::field::FieldElement51;
use crate::backend::vector::avx2::constants::{
    P_TIMES_2_HI, P_TIMES_2_LO,
};

use curve25519_dalek_derive::unsafe_target_feature;





/// The `Lanes` enum represents a subset of the lanes `A,B,C,D` of a
/// `FieldElement2625x4`.
///
/// It's used to specify blend operations without
/// having to know details about the data layout of the
/// `FieldElement2625x4`.
#[allow(clippy::upper_case_acronyms)]
#[derive(Copy, Clone, Debug)]
pub(crate) enum Lanes {
    C,
    D,
    AB,
    AC,
    CD,
    AD,
    BC,
    ABCD,
}

/// The `Shuffle` enum represents a shuffle of a `FieldElement2625x4`.
///
/// The enum variants are named by what they do to a vector \\(
/// (A,B,C,D) \\); for instance, `Shuffle::BADC` turns \\( (A, B, C,
/// D) \\) into \\( (B, A, D, C) \\).
#[allow(clippy::upper_case_acronyms)]
#[derive(Copy, Clone, Debug)]
pub(crate) enum Shuffle {
    AAAA,
    BBBB,
    CACA,
    DBBD,
    ADDA,
    CBCB,
    ABAB,
    BADC,
    BACD,
    ABDC,
}

/// A vector of four field elements.
///
/// Each operation on a `FieldElement2625x4` has documented effects on
/// the bounds of the coefficients.  This API is designed for speed
/// and not safety; it is the caller's responsibility to ensure that
/// the post-conditions of one operation are compatible with the
/// pre-conditions of the next.
#[derive(Clone, Copy, Debug)]
pub(crate) struct FieldElement2625x4(pub(crate) [u32x8; 5]);


#[unsafe_target_feature("avx2")]
impl FieldElement2625x4 {
    pub(crate) const ZERO: FieldElement2625x4 = FieldElement2625x4([u32x8::splat_const::<0>(); 5]);

    /// Split this vector into an array of four (serial) field
    /// elements.
    #[rustfmt::skip] // keep alignment of extracted lanes
    pub(crate) fn split(&self) -> [FieldElement51; 4] {
        let mut out = [FieldElement51::ZERO; 4];
        for i in 0..5 {
            let a_2i   = self.0[i].extract::<0>() as u64; //
            let b_2i   = self.0[i].extract::<1>() as u64; //
            let a_2i_1 = self.0[i].extract::<2>() as u64; // `.
            let b_2i_1 = self.0[i].extract::<3>() as u64; //  | pre-swapped to avoid
            let c_2i   = self.0[i].extract::<4>() as u64; //  | a cross lane shuffle
            let d_2i   = self.0[i].extract::<5>() as u64; // .'
            let c_2i_1 = self.0[i].extract::<6>() as u64; //
            let d_2i_1 = self.0[i].extract::<7>() as u64; //

            out[0].0[i] = a_2i + (a_2i_1 << 26);
            out[1].0[i] = b_2i + (b_2i_1 << 26);
            out[2].0[i] = c_2i + (c_2i_1 << 26);
            out[3].0[i] = d_2i + (d_2i_1 << 26);
        }

        out
    }

    /// Given \\((A,B,C,D)\\), compute \\((-A,-B,-C,-D)\\), without
    /// performing a reduction.
    ///
    /// # Preconditions
    ///
    /// The coefficients of `self` must be bounded with \\( b < 0.999 \\).
    ///
    /// # Postconditions
    ///
    /// The coefficients of the result are bounded with \\( b < 1 \\).
    #[inline]
    pub(crate) fn negate_lazy(&self) -> FieldElement2625x4 {
        // The limbs of self are bounded with b < 0.999, while the
        // smallest limb of 2*p is 67108845 > 2^{26+0.9999}, so
        // underflows are not possible.
        FieldElement2625x4([
            P_TIMES_2_LO - self.0[0],
            P_TIMES_2_HI - self.0[1],
            P_TIMES_2_HI - self.0[2],
            P_TIMES_2_HI - self.0[3],
            P_TIMES_2_HI - self.0[4],
        ])
    }
}
