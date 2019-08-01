import re
import copy
from collections import OrderedDict


GT_LINE = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
STATS_HEADER = ['##INFO=<ID=Identified,Number=1,Type=String,Description="Source \
VCF for the merged record">\n',
                '##INFO=<ID=AD_mean,Number=1,Type=Float,Description="Mean of \
allelic depths in source vcfs">\n',
                '##INFO=<ID=AD_sd,Number=1,Type=Float,Description="Standard \
deviation of allelic depths in source vcfs">\n',
                '##INFO=<ID=DP_mean,Number=1,Type=Float,Description="Mean of \
read depths in source vcfs">\n',
                '##INFO=<ID=DP_sd,Number=1,Type=Float,Description="Standard \
deviation of read depths in source vcfs">\n',
                '##FORMAT=<ID=AD,Number=1,Type=Float,Description="Mean of \
allelic depths in source vcfs">\n',
                '##FORMAT=<ID=DP,Number=1,Type=Float,Description="Mean of read \
depths in source vcfs">\n']
SOMATIC_STATS_HEADER = ['##INFO=<ID=Identified,Number=1,Type=String, \
Description="Source VCF for the merged record">\n',
                        '##INFO=<ID=AD_mean_normal,Number=1,Type=Float, \
Description="Mean of allelic depths in source vcfs">\n',
                        '##INFO=<ID=AD_sd_normal,Number=1,Type=Float, \
Description="Standard deviation of allelic depths in source vcfs">\n',
                        '##INFO=<ID=AD_mean_tumor,Number=1,Type=Float, \
Description="Mean of allelic depths in source vcfs">\n',
                        '##INFO=<ID=AD_sd_tumor,Number=1,Type=Float, \
Description="Standard deviation of allelic depths in source vcfs">\n',
                        '##INFO=<ID=DP_mean_normal,Number=1,Type=Float, \
Description="Mean of read depths in source vcfs">\n',
                        '##INFO=<ID=DP_sd_normal,Number=1,Type=Float, \
Description="Standard deviation of read depths in source vcfs">\n',
                        '##INFO=<ID=DP_mean_tumor,Number=1,Type=Float, \
Description="Mean of read depths in source vcfs">\n',
                        '##INFO=<ID=DP_sd_tumor,Number=1,Type=Float, \
Description="Standard deviation of read depths in source vcfs">\n',
                        '##FORMAT=<ID=AD,Number=1,Type=Float, \
Description="Mean of allelic depths in source vcfs">\n',
                        '##FORMAT=<ID=DP,Number=1,Type=Float, \
Description="Mean of read depths in source vcfs">\n']
HEADER = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'

#####################################################################


def extract_cols(i_dict, f_dict, cols):
    """Select the columns from INFO/FORMAT using list 'cols'
    Input:
        info and format dictionaries of header object
        {column_name: header object}
        cols, a list of column names specified by the user
    Output:
        new info and format dictionaries
    """
    i_cols, f_cols = OrderedDict(), OrderedDict()
    for col in cols:
        i_header, f_header = i_dict.get(col), f_dict.get(col)
        if i_header:
            i_cols.update({col: i_header})
        #if f_header and (f_header.meta_id in ['AD', 'DP']):
        # Move all the FORMAT columns (except GT) to INFO
        if f_header and (f_header.meta_id != "GT"):
            f_header.meta = 'INFO'
            i_cols.update({col: f_header})
        # If "GT" specified, keep the GT in the FORMAT
        elif f_header and (f_header.meta_id == "GT"):
            f_cols.update({col: f_header})
    return i_cols, f_cols


def extract_cols_somatic(i_dict, f_dict, cols):
    """Select the columns from INFO/FORMAT using list 'cols'
    Input:
        info and format dictionaries of header object
        {column_name: header object}
        cols, a list of column names specified by the user
    Output:
        new info and format dictionaries
    """

    i_cols, f_cols = OrderedDict(), OrderedDict()
    cols_ = [i for i in cols if i != 'GT']
    for col in cols:
        i_header, f_header = i_dict.get(col), f_dict.get(col)
        if i_header:
            i_cols.update({col: i_header})
        # Put the caller's AD and DP in INFO
        if f_header and f_header.meta_id in cols_:
            f_header.meta = 'INFO'
            # If this column is coming from FORMAT, it would have normal/tumor
            # value
            normal_header = f_header
            tumor_header = copy.deepcopy(f_header)
            normal_header.meta_id += '_normal'
            tumor_header.meta_id += '_tumor'
            i_cols.update({col + '_normal': normal_header})
            i_cols.update({col + '_tumor': tumor_header})

        elif f_header:
            f_cols.update({col: f_header})
    # Force all vcf to have GT
    if 'GT' in cols and ('GT' not in f_cols.keys()):
        f_cols.update({'GT': VcfHeader(GT_LINE)})
    return i_cols, f_cols

#####################################################################


class VcfHeader:
    """Meta info lines in vcf
    """

    def __init__(self, line):
        """Create object from the meta informaiton line
        """

        line = re.split('<|>|"', line)
        h_dict = {}
        for i in re.split(',(?!\s)', line[1]):
            h_dict.update({i.split('=')[0]: i.split('=')[1]})

        self.meta = line[0][:-1].replace('#', '')
        self.meta_id = h_dict.get('ID', None)
        self.meta_number = h_dict.get('Number', None)
        self.meta_type = h_dict.get('Type', None)
        self.meta_description = line[2]
        self.normal_tumor = False

    def add_caller(self, caller):
        """Input: header object
           Output: add the caller to the header. If comming from somatic
           FORMAT, will create normal and tumor header line at the same time.
        """

        if self.meta_id not in ['AC', 'AN']:
            self.meta_id += '_{}'.format(caller)
            self.meta_description = (self.meta_description +
                                     ' ({})'.format(caller))
        return self

    def write(self):
        """Write modified meta informaiton line
        """

        return ('##{}=<ID={},Number={},Type={},Description="{}">\n'
                .format(self.meta, self.meta_id, self.meta_number,
                        self.meta_type, self.meta_description))

