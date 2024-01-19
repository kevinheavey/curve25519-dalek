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

//! This module contains constants used by the AVX2 backend.

use crate::backend::vector::packed_simd::u32x8;

/// The low limbs of (2p, 2p, 2p, 2p), so that
/// ```ascii,no_run
/// (2p, 2p, 2p, 2p) = [P_TIMES_2_LO, P_TIMES_2_HI, P_TIMES_2_HI, P_TIMES_2_HI, P_TIMES_2_HI]
/// ```
pub(crate) static P_TIMES_2_LO: u32x8 = u32x8::new_const(
    67108845 << 1,
    67108845 << 1,
    33554431 << 1,
    33554431 << 1,
    67108845 << 1,
    67108845 << 1,
    33554431 << 1,
    33554431 << 1,
);

/// The high limbs of (2p, 2p, 2p, 2p), so that
/// ```ascii,no_run
/// (2p, 2p, 2p, 2p) = [P_TIMES_2_LO, P_TIMES_2_HI, P_TIMES_2_HI, P_TIMES_2_HI, P_TIMES_2_HI]
/// ```
pub(crate) static P_TIMES_2_HI: u32x8 = u32x8::new_const(
    67108863 << 1,
    67108863 << 1,
    33554431 << 1,
    33554431 << 1,
    67108863 << 1,
    67108863 << 1,
    33554431 << 1,
    33554431 << 1,
);

/// The low limbs of (16p, 16p, 16p, 16p), so that
/// ```ascii,no_run
/// (16p, 16p, 16p, 16p) = [P_TIMES_16_LO, P_TIMES_16_HI, P_TIMES_16_HI, P_TIMES_16_HI, P_TIMES_16_HI]
/// ```
pub(crate) static P_TIMES_16_LO: u32x8 = u32x8::new_const(
    67108845 << 4,
    67108845 << 4,
    33554431 << 4,
    33554431 << 4,
    67108845 << 4,
    67108845 << 4,
    33554431 << 4,
    33554431 << 4,
);

/// The high limbs of (16p, 16p, 16p, 16p), so that
/// ```ascii,no_run
/// (16p, 16p, 16p, 16p) = [P_TIMES_16_LO, P_TIMES_16_HI, P_TIMES_16_HI, P_TIMES_16_HI, P_TIMES_16_HI]
/// ```
pub(crate) static P_TIMES_16_HI: u32x8 = u32x8::new_const(
    67108863 << 4,
    67108863 << 4,
    33554431 << 4,
    33554431 << 4,
    67108863 << 4,
    67108863 << 4,
    33554431 << 4,
    33554431 << 4,
);
