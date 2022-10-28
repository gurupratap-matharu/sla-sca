# Snow Line Altitude Monitoring (Argentina)

This project processes geospatial data of glaciers in Argentina and identifies the
behaviour of level of snow line altitude across a time period.

## Features

- Uses Earth Engine data catalog
- Is generic and can work for any glacier in the world
- Uses 100% python for processing earth engine data
- Can act as a server for a web application

## Contributing

Contributions are always welcome!

See `contributing.md` for ways to get started.

Please adhere to this project's `code of conduct`.

## Lessons Learned

- How does Earth Engine work?
- What is raster data, Image Collections, Feature Collections?
- How to run machine learning model on Geospatial datasets?
- How to use Python API for earth engine?
- How to visualize Earth engine data on a Map with layers?

## Run Locally

Clone the project

```bash
  git clone https://github.com/gurupratap-matharu/sla-sca
```

Go to the project directory

```bash
  cd sla-sca
```

Install dependencies

```bash
  poetry install
```

Process the data

```bash
  python main.py
```

## Authors

- [@gurupratap-matharu](https://www.github.com/gurupratap-matharu)
