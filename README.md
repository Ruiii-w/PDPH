# PDPH

## Env

- gflags
- GCC

## Structure

**Header files:**

`PDPH.h` is the header file of our proposal for testing.

`container.h` contains definitions for various basic index and storage structures.

`murmur3.h` is the header file for the hash function used in PDPH, with MurmurHash3 as the default.

`util.h` includes numerous macro definitions and helper functions.


## usage

```
./PDPH -index=[index name] -dataset=[dataset file path] -thread=[number of threads]
```

