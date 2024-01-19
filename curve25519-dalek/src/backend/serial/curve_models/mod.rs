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

//! Internal curve representations which are not part of the public API.
//!
//! # Curve representations
//!
//! Internally, we use several different models for the curve.  Here
//! is a sketch of the relationship between the models, following [a
//! post][smith-moderncrypto]
//! by Ben Smith on the `moderncrypto` mailing list.  This is also briefly
//! discussed in section 2.5 of [_Montgomery curves and their
//! arithmetic_][costello-smith-2017] by Costello and Smith.
//!
//! Begin with the affine equation for the curve,
//! $$
//!     -x\^2 + y\^2 = 1 + dx\^2y\^2.
//! $$
//! Next, pass to the projective closure \\(\mathbb P\^1 \times \mathbb
//! P\^1 \\) by setting \\(x=X/Z\\), \\(y=Y/T.\\)  Clearing denominators
//! gives the model
//! $$
//!     -X\^2T\^2 + Y\^2Z\^2 = Z\^2T\^2 + dX\^2Y\^2.
//! $$
//! In `curve25519-dalek`, this is represented as the `CompletedPoint`
//! struct.
//! To map from \\(\mathbb P\^1 \times \mathbb P\^1 \\), a product of
//! two lines, to \\(\mathbb P\^3\\), we use the [Segre
//! embedding](https://en.wikipedia.org/wiki/Segre_embedding)
//! $$
//!     \sigma : ((X:Z),(Y:T)) \mapsto (XY:XT:ZY:ZT).
//! $$
//! Using coordinates \\( (W_0:W_1:W_2:W_3) \\) for \\(\mathbb P\^3\\),
//! the image \\(\sigma (\mathbb P\^1 \times \mathbb P\^1) \\) is the
//! surface defined by \\( W_0 W_3 = W_1 W_2 \\), and under \\(
//! \sigma\\), the equation above becomes
//! $$
//!     -W\_1\^2 + W\_2\^2 = W\_3\^2 + dW\_0\^2,
//! $$
//! so that the curve is given by the pair of equations
//! $$
//! \begin{aligned}
//!     -W\_1\^2 + W\_2\^2 &= W\_3\^2 + dW\_0\^2, \\\\  W_0 W_3 &= W_1 W_2.
//! \end{aligned}
//! $$
//! Up to variable naming, this is exactly the "extended" curve model
//! introduced in [_Twisted Edwards Curves
//! Revisited_][hisil-wong-carter-dawson-2008] by Hisil, Wong, Carter,
//! and Dawson.  In `curve25519-dalek`, it is represented as the
//! `EdwardsPoint` struct.  We can map from \\(\mathbb P\^3 \\) to
//! \\(\mathbb P\^2 \\) by sending \\( (W\_0:W\_1:W\_2:W\_3) \\) to \\(
//! (W\_1:W\_2:W\_3) \\).  Notice that
//! $$
//!     \frac {W\_1} {W\_3} = \frac {XT} {ZT} = \frac X Z = x,
//! $$
//! and
//! $$
//!     \frac {W\_2} {W\_3} = \frac {YZ} {ZT} = \frac Y T = y,
//! $$
//! so this is the same as if we had started with the affine model
//! and passed to \\( \mathbb P\^2 \\) by setting \\( x = W\_1 / W\_3
//! \\), \\(y = W\_2 / W\_3 \\).
//! Up to variable naming, this is the projective representation
//! introduced in in [_Twisted Edwards
//! Curves_][bernstein-birkner-joye-lange-peters-2008] by Bernstein,
//! Birkner, Joye, Lange, and Peters.  In `curve25519-dalek`, it is
//! represented by the `ProjectivePoint` struct.
//!
//! # Passing between curve models
//!
//! Although the \\( \mathbb P\^3 \\) model provides faster addition
//! formulas, the \\( \mathbb P\^2 \\) model provides faster doubling
//! formulas.  Hisil, Wong, Carter, and Dawson therefore suggest mixing
//! coordinate systems for scalar multiplication, attributing the idea
//! to [a 1998 paper][cohen-miyaji-ono-1998] of Cohen, Miyagi, and Ono.
//!
//! Their suggestion is to vary the formulas used by context, using a
//! \\( \mathbb P\^2 \rightarrow \mathbb P\^2 \\) doubling formula when
//! a doubling is followed
//! by another doubling, a \\( \mathbb P\^2 \rightarrow \mathbb P\^3 \\)
//! doubling formula when a doubling is followed by an addition, and
//! computing point additions using a \\( \mathbb P\^3 \times \mathbb P\^3
//! \rightarrow \mathbb P\^2 \\) formula.
//!
//! The `ref10` reference implementation of [Ed25519][ed25519], by
//! Bernstein, Duif, Lange, Schwabe, and Yang, tweaks
//! this strategy, factoring the addition formulas through the
//! completion \\( \mathbb P\^1 \times \mathbb P\^1 \\), so that the
//! output of an addition or doubling always lies in \\( \mathbb P\^1 \times
//! \mathbb P\^1\\), and the choice of which formula to use is replaced
//! by a choice of whether to convert the result to \\( \mathbb P\^2 \\)
//! or \\(\mathbb P\^3 \\).  However, this tweak is not described in
//! their paper, only in their software.
//!
//! Our naming for the `CompletedPoint` (\\(\mathbb P\^1 \times \mathbb
//! P\^1 \\)), `ProjectivePoint` (\\(\mathbb P\^2 \\)), and
//! `EdwardsPoint` (\\(\mathbb P\^3 \\)) structs follows the naming in
//! Adam Langley's [Golang ed25519][agl-ed25519] implementation, which
//! `curve25519-dalek` was originally derived from.
//!
//! Finally, to accelerate readditions, we use two cached point formats
//! in "Niels coordinates", named for Niels Duif,
//! one for the affine model and one for the \\( \mathbb P\^3 \\) model:
//!
//! * `AffineNielsPoint`: \\( (y+x, y-x, 2dxy) \\)
//! * `ProjectiveNielsPoint`: \\( (Y+X, Y-X, Z, 2dXY) \\)
//!
//! [smith-moderncrypto]: https://moderncrypto.org/mail-archive/curves/2016/000807.html
//! [costello-smith-2017]: https://eprint.iacr.org/2017/212
//! [hisil-wong-carter-dawson-2008]: https://www.iacr.org/archive/asiacrypt2008/53500329/53500329.pdf
//! [bernstein-birkner-joye-lange-peters-2008]: https://eprint.iacr.org/2008/013
//! [cohen-miyaji-ono-1998]: https://link.springer.com/content/pdf/10.1007%2F3-540-49649-1_6.pdf
//! [ed25519]: https://eprint.iacr.org/2011/368
//! [agl-ed25519]: https://github.com/agl/ed25519

#![allow(non_snake_case)]

use crate::field::FieldElement;

/// A pre-computed point in the affine model for the curve, represented as
/// \\((y+x, y-x, 2dxy)\\) in "Niels coordinates".
///
/// More details on the relationships between the different curve models
/// can be found in the module-level documentation.
// Safe to derive Eq because affine coordinates.
#[derive(Copy, Clone, Eq, PartialEq)]
#[allow(missing_docs)]
pub(crate) struct AffineNielsPoint {
    pub y_plus_x: FieldElement,
    pub y_minus_x: FieldElement,
    pub xy2d: FieldElement,
}
