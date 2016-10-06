# pypie
A suite of functions to calibrate mass spectral data and extract photo ionization efficiency (PIE) curves from their energy series. This script was designed for use with data collected at the Chemical Dynamics 9.0.2 beam line at the Berkeley Labs Advanced Light Source but may be used with any mass spectral series data after inheriting from the Pie class and redefining the load method. Pypie requires the numpy, matplotlib, and uncertainties libraries. It is also recommended to install the periodictable library for use of the formula function but this is not required. See the example uses below. 

## Example Usage
The suite may be tested by downloading the included test data set "premix_data.zip". It contains 3 repeat scans for demonstration of both the Pie and uPie (PIE with uncertainty) classes. For Pie, choose any data set. For uPie, use all three. 

### The Pie object
```python
    from pypie import Pie
    from periodictable import formula # If not using formula, masses should be specified explicitly.
    
    path_to_data = '/path/to/data/premix_1.dat'
    pie = Pie()
    pie.load(path_to_data)
    m1, m2 = formula('He').mass, formula('C6H6').mass
    t1, t2 = 250, 8000
```
