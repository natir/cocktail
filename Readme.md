# Cocktail

[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/natir/cocktail/blob/master/LICENSE)
![CI](https://github.com/natir/cocktail/workflows/CI/badge.svg)
[![Documentation](https://github.com/natir/cocktail/workflows/Documentation/badge.svg)](https://natir.github.io/cocktail/cocktail)
[![CodeCov](https://codecov.io/gh/natir/cocktail/branch/master/graph/badge.svg)](https://codecov.io/gh/natir/cocktail)

Cocktail it's a rust crate, python module, c library, to convert DNA in kmer 2 bit representation and get is cannonical version.

**Warning this isn't stable, API can change any time**

- [Usage](#usage)
- [Minimum supported Rust version](#minimum-supported-rust-version)
- [Citation](#citation)

# Usage

## Rust

In `[dependencies]` section of your `Cargo.toml` add this: 
```
cocktail = { git="https://github.com/natir/cocktail.git" }
```

## Python 

You need [rust toolchain setup on your system](https://rustup.rs/)

Give this to pip:
```
git+https://github.com/natir/cocktail.git#egg=cocktail&subdirectory=dist/python
```

## C

You need [rust toolchain setup on your system](https://rustup.rs/)

```
git clone https://github.com/natir/cocktail.git
cd cocktail
cargo build --release

cbindgen --config cbindgen.toml --crate cocktail --output dist/c/cocktail.h
cd dist/c/
make
./test
```

Dynamic and static library is avaible her `target/release/libcocktail.{a|so}` header is her `dist/c/cocktail.h`. To build a C programe you need to add `-lpthread -lm -ldl` durring linking phase.

## Minimum supported Rust version

Currently the minimum supported Rust version is 1.45.0.

## Citation

WIP
