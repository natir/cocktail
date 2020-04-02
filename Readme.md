# Cocktail

Cocktail it's a rust crate, python module, c library, to convert DNA in kmer 2 bit representation and get is cannonical version.

**Warning this isn't stable API can change any time**

# Usage

## Rust

In `[dependencies]` section of your `Cargo.toml` add this: 
```
cocktail = { git="https://github.com/natir/cocktail.git" }
```

## Python 

You need [rust toolchain setup on your system](https://rustup.rs/)

In `requirements.txt` add this:
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
