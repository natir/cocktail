all: test.c ../../target/release/libcocktail.a
	gcc -Wall test.c ../../target/release/libcocktail.a -lpthread -lm -ldl -o test

../../target/release/libcocktail.a:
	cd ../..
	cargo build --release
