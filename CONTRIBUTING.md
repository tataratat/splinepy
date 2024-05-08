# Contributing
splinepy welcomes and appreciates discussions, issues and pull requests!

## Quick start
Once the repo is forked, one possible starting point would be creating a new python environments, for example, using [conda](https://docs.conda.io/en/latest/miniconda.html) with, for example, `python=3.10`
```bash
conda create -n splinepyenv python=3.10
conda activate splinepyenv
git clone git@github.com:<path-to-your-fork>
cd splinepy  # or <forkname>
git submodule update --init --recursive
git checkout -b new-feature0
pip install -e. -v --config-settings=cmake.args=-DSPLINEPY_MORE=OFF --config-settings=cmake.build-type="Debug"
```
`--config-settings=cmake.args=-DSPLINEPY_MORE` build argument builds splines up to 3D (both parametric and physical dimensions if they are part of template parameters), and that way we can reduce compile time. `--config-settings=cmake.build-type="Debug"` build also reduces compile time. We are experimenting with ways to reduce compile time during development. Let us know if you have a good idea!

## Python style / implementation preferences
- use [`PEP 8`](https://peps.python.org/pep-0008/) style guide
- no complex comphrehensions: preferably fits in a line, 2 lines max if it is totally necessary
- use first letter abbreviations in element loops:  `for kv in knot_vectors`
- use `i`, `j`, `k`, `l` for pure index: `for i, kv in enumerate(knot_vectors)`
- module local import with a leading underscore: `from splinepy import settings as _settings`
- make sure function wrappers are properly documented (see, e.g., helpme.create.Creator's member functions)

## C++ style / implementation preferences
For c++, we've prepared a `.clang-format`, with that you can just run `clang-format`. We closely follow naming scheme suggested by [google stype guide](https://google.github.io/styleguide/cppguide.html#Naming), with a clear exception of file naming: we use [`.hpp`, `.inl`, `.cpp`].
Here's another preference:
- `#pragma once`

## Automatic formatting / style check
To check the format and style of your code use the following commands:
```bash
pip install pre-commit
precommit run -a
```

## Local documentation build
```bash
pip install -r docs/requirements.txt
python3 docs/source/handle_markdown.py
python3 docs/markdown/extra_docs.py
sphinx-build -b html docs/source docs/build -E
# you can now open `docs/build/index.html` with your browser
```
