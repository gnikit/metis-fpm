# METIS for fpm

This provides a Fortran API and a source repackaging of the
[METIS](https://github.com/KarypisLab/METIS) library originally developed by Karypis Lab

## TODO

- [ ] Run tests in CI
- [ ] Make CI open a PR for changes if the repository is dirty
- [ ] Add Fortran API procedural and OOP interfaces

## Usage

```sh
fpm build --c-flag "-DIDXTYPEWIDTH=64 -DREALTYPEWIDTH=64"
```

To use `metis` as a dependency in your `fpm` project, add the following to your `fpm.toml` file:

```toml
[dependencies]
gklib = { git = "https://github.com/gnikit/metis-fpm.git" }
```

## License
