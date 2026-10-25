"""Read centroided LC-MS spectra."""

from pyteomics import mzml, mzxml


class MS1Chromatogram:
    """
    A class to represent MS1 Chromatograms.

    Attributes:
    - rt_list (list): List of retention times.
    - spectra (list): List of spectra.
    """

    def __init__(self):
        self.rt_list = []
        self.spectra = []


def load_rawdata(rawdata_file):
    ms_chrom = MS1Chromatogram()
    if rawdata_file.endswith(".mzXML"):
        rawdata_reader = mzxml.MzXML(rawdata_file)
    else:
        rawdata_reader = mzml.MzML(rawdata_file)
    for spectrum in rawdata_reader:
        ms_chrom.rt_list.append(spectrum["retentionTime"] * 60)
        intensity_array = spectrum["intensity array"]
        mz_array = spectrum["m/z array"][intensity_array > 0]
        intensity_array = intensity_array[intensity_array > 0]
        ms_chrom.spectra.append((mz_array, intensity_array))

    return ms_chrom
