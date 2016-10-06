# pypie
A suite of functions to calibrate mass spectral data and extract photo ionization efficiency (PIE) curves from their energy series. This script was designed for use with data collected at the Chemical Dynamics 9.0.2 beam line at the Berkeley Labs Advanced Light Source but may be used with any mass spectral series data after inheriting from the Pie class and redefining the load method. Pypie requires the numpy, matplotlib, and uncertainties libraries. 

## Example Usage

    from pypie import Pie
    paths_to_data = (
        '/path/to/data/repeat_1',
        '/path/to/data/repeat_2',
        '/path/to/data/repeat_3',
        ) # These paths should be repeats for 
    pie = Pie(paths_to_data)
     
