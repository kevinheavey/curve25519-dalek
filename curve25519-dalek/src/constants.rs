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
//! Various constants, such as the Ristretto and Ed25519 basepoints.

#![allow(non_snake_case)]

use cfg_if::cfg_if;



cfg_if! {
    if #[cfg(curve25519_dalek_backend = "fiat")] {
        #[cfg(curve25519_dalek_bits = "32")]
        pub(crate) use crate::backend::serial::fiat_u32::constants::*;
        #[cfg(curve25519_dalek_bits = "64")]
        pub(crate) use crate::backend::serial::fiat_u64::constants::*;
    } else {
        #[cfg(curve25519_dalek_bits = "32")]
        pub(crate) use crate::backend::serial::u32::constants::*;
        #[cfg(curve25519_dalek_bits = "64")]
        pub(crate) use crate::backend::serial::u64::constants::*;
    }
}
