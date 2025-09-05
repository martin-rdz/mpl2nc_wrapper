
# mpl2nc wrapper

Simple wrapper for [mpl2nc](https://github.com/peterkuma/mpl2nc) that produces daily netcdf files which also include backscatter and depolarization.

```
python3 wrapper_mpl2nc_cloudnet.py --date 20240101
```

The calibration constant can be estimated in (not too low) liquid clouds [O'Conner et al., 2004](https://journals.ametsoc.org/view/journals/atot/21/5/1520-0426_2004_021_0777_atfaoc_2_0_co_2.xml), but is currently hard coded.
Deadtime correction has to be taken from the table in the manual as a csv file. The binary file given might have numerical inaccuracies.

The resulting daily netcdf files are quite large (~2.8GB), but might contain superfluous information (see [output_nc_dump.txt](output_nc_dump.txt)).

## Custom overlap functions

A customly derived overlap function can also be stored in the binary format needed for the library:

```
import struct

def read_overlap(filename):
    """from mpl2nc https://github.com/peterkuma/mpl2nc/blob/94c92358c939ebb4ba13a71674f4a18695de051a/mpl2nc.py#L208"""

    with open(filename, 'rb') as f:
        buf = f.read()
        n = len(buf)//16
        f.seek(0)
        d = {'ol_number_bins': np.array(n, np.uint32)}
        for x in ['ol_range', 'ol_overlap']:
            buf = f.read(8*n)
            if len(buf) < 8*n:
                raise IOError('incomplete %s data' % x)
            a = struct.unpack_from('<' + 'd'*n, buf)
            d[x] = np.array(a, np.float64)
    return d

manufacturer_overlap = read_overlap('MMPL5081_Overlap_202203111400.bin')
man_range = manufacturer_overlap['ol_range']
```


```
def write_overlap(filename, d):
    """Write binary data in the same format as read_overlap"""
    ol_range = d['ol_range']
    ol_overlap = d['ol_overlap']
    if len(ol_range) != len(ol_overlap):
        raise ValueError("ol_range and ol_overlap must have the same length")
    n = len(ol_range)
    with open(filename, 'wb') as f:
        # Write ol_range
        f.write(struct.pack('<' + 'd'*n, *ol_range))
        # Write ol_overlap
        f.write(struct.pack('<' + 'd'*n, *ol_overlap))

write_overlap('MMPL5081_Overlap_202203111400_outputtest.bin', manufacturer_overlap)
```