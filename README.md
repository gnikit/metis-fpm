# METIS for fpm

This provides a Fortran API and a source repackaging of the
[METIS](https://github.com/KarypisLab/METIS) library originally developed by Karypis Lab

## TODO

- [ ] Run tests in CI
- [ ] Make CI open a PR for changes if the repository is dirty
- [ ] Add Fortran API procedural and OOP interfaces

## Usage

```sh
fpm build --flag "-DIDXTYPEWIDTH=64 -DREALTYPEWIDTH=64" --c-flag "-DIDXTYPEWIDTH=64 -DREALTYPEWIDTH=64"
```

To use `metis` as a dependency in your `fpm` project, add the following to your `fpm.toml` file:

```toml
[dependencies]
gklib = { git = "https://github.com/gnikit/metis-fpm.git" }
```

## License

> MIT License

[metis-fpm](https://github.com/gnikit/metis-fpm) is distributed under the MIT License. See [LICENSE](LICENSE) for more information.
[METIS](https://github.com/KarypisLab/METIS) is distributed under the terms of the Apache License, Version 2.0
