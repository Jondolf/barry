# Barry

> Note: The name Barry is not final and is subject to change.

**Barry** is a fork of **[Parry]**, providing a set of 2 and 3-dimensional geometric and collision detection libraries
using [Glam] math types through [`bevy_math`].
These libraries are `barry2d`, `barry3d`, `barry2d-f64`, and `barry3d-f64`. It is forever free
and open-source!

## Warning

Barry is currently at a very early stage, and the migration to Bevy's math types is still on-going.
It is not usable yet, and relies on the main branches of some dependencies or even local changes that I have not
released yet. I will update this once the core migration is fully done and the functionality of the original Parry
is working as intended.

## Motivation

The Rust ecosystem does not have a good collision detection library that uses [Glam]. The most popular library,
[Parry], uses [Nalgebra].

Nalgebra is an incredible linear algebra library that is a good fit for things like graphics and physics simulations,
but many major projects, namely the [Bevy game engine][Bevy], use Glam due to reasons like its ease of use and minimal nature.
Having two different math libraries in the same project can be confusing and inconvenient for both users and contributors,
and it also increases compile times and dependency size. Ideally, there would only be need for a single math library.

This is where Barry (name not final) comes in! It aims to provide the Glam-based side of the ecosystem a fully featured
and mature collision detection library. Instead of reinventing the wheel and creating everything from scratch,
Barry is a fork of the brilliant collision detection crate Parry, modified to use Glam types.

Specifically, Barry uses [`bevy_math`], which uses Glam. A core goal of Barry is to design the ideal collision detection
crate for the Bevy game engine that aligns with its needs and philosophies. For good interoperability, Barry will use
Bevy's geometric primitives (currently on the main branch, not released) instead of defining its own separate set of shapes.

Note that this does *not* mean that you need to use Bevy itself; Barry simply uses Bevy's math types, and can be used fully
independently of Bevy.

## Goals

Below are some of Barry's long-term goals.

- Provide a Glam-based alternative to Parry (mostly done)
- Develop a collision detection library with great Bevy interoperability
  - Use Glam
  - Use Bevy's geometric primitives for collider shapes
  - Use Bevy's bounding volumes
  - Use Bevy's ray types
  - Support using Bevy's threadpool instead of [`rayon`] to unblock threads when using Bevy
- Improve the documentation and examples, hopefully contributing back to Parry
- Have freedom to add new features and improvements more efficiently, hopefully contributing back to Parry
- (Personal goal) Get a better understanding of the details of collision detection

[Parry]: https://parry.rs/
[Nalgebra]: https://nalgebra.org/
[Glam]: https://github.com/bitshifter/glam-rs
[`bevy_math`]: https://docs.rs/bevy_math/latest/bevy_math/
[Bevy]: https://bevyengine.org/
[`rayon`]: https://github.com/rayon-rs/rayon

## License

Barry is free, open source and permissively licensed!
Except where noted (below and/or in individual files), all code in this repository is licensed under the
Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)).

Some of the engine's code carries additional copyright notices and license terms due to their external origins.
These files will contain a note of the separate licensing, typically at the top of the file or next to the code where it is used.

### Your contributions

Unless you explicitly state otherwise,
any contribution intentionally submitted for inclusion in the work by you,
as defined in the Apache-2.0 license,
shall be licensed as above,
without any additional terms or conditions.
