import sys
from collections import OrderedDict as od
from vcfheader import VcfHeader, extract_cols, extract_cols_somatic
from variant import Variant

##############################################################################

AF_LINE = '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, \
for each ALT allele, in the same order as listed">'

##############################################################################


def sort_vcf(priority, callers, vcf_list):
    """Sort the vcf
    Output two list with the matched sorting order:
        vcf_list
        callers
    """
    if priority:
        # Sort in user specified order
        if set(priority) != set(callers):
            sys.exit("The priority '{}' does not match the callers list '{}' \
from the vcfs".format(' '.join(priority), ' '.join(callers)))
        else:
            priority_dict = {}
            for i, caller in enumerate(priority):
                priority_dict[caller] = i
            vcf_list = sorted(vcf_list, key=lambda vcf:
                              priority_dict[vcf.caller])
            callers = sorted(callers, key=lambda caller: priority_dict[caller])
    else:
        # sort in alphatic order
        vcf_list = sorted(vcf_list, key=lambda vcf: vcf.caller.lower())
        callers = sorted(callers, key=lambda caller: caller.lower())
    return vcf_list, callers

##############################################################################


class NormalisedVcf:
    """Extract and normalise selected columns from a vcf
    Attributes:
    name: string, name of the vcf file
    caller: string, extracted from the vcf line ##source
    meta_info: meta info lines (list)
    header: #CHROM\t.... (string)
    variants: variant objects query by 'CHROM_POS_ID_REF_ALF' (dictionary)
    """

    def __init__(self, vcfname):
        self.name = vcfname
        self.caller = ''
        self.meta_info = []
        self.header = ''
        self.variants = {}

    def process_vcf(self, cols):
        """Build object from vcf
        """
        vcf = open(self.name, 'r')
        info_dict, format_dict = {}, {}

        # Read the meta-information lines from the vcf
        for i, line in enumerate(vcf):
            # Handle exceptions: the AF will be calcualted regardless;
            if line.startswith('##FORMAT=<ID=AF'):
                pass
            # Select the INFO/FORMAT lines
            elif line.startswith('##FORMAT'):
                vcf_header = VcfHeader(line)
                format_dict.update({vcf_header.meta_id: vcf_header})
            elif line.startswith('##INFO'):
                vcf_header = VcfHeader(line)
                info_dict.update({vcf_header.meta_id: vcf_header})
            # Keep other meta-info lines
            elif line.startswith('##'):
                if line.startswith('##source='):
                    self.caller = line.replace('##source=', '').strip()
                self.meta_info.append(line)
            else:
                break
            # Only extract the (filtered) DP in the format
            if "DP" in info_dict.keys() and format_dict.keys():
                info_dict.pop("DP", None)

        if not self.caller:
            sys.exit("Cannot identify caller from file {}\nPlease add caller \
                     identify line '##source=(caller name)' to vcf header"
                     .format(self.name))

        # When user specify the AF and vcf does not have, try to calculate that 
        # for the user
        if ('AF' in cols) and ('AF' not in info_dict.keys()):
            vcf_header = VcfHeader(AF_LINE)
            info_dict.update({vcf_header.meta_id: vcf_header})

        # Select the columns from INFO/FORMAT
        info_cols, format_cols = extract_cols(info_dict, format_dict, cols)

        # Add the INFO line (with caller) / FORMAT (unchanged) to header_list
        self.meta_info += [VcfHeader.write(VcfHeader.add_caller(v,
                           self.caller)) for k, v in info_cols.items()]
        self.meta_info += [VcfHeader.write(v) for k, v in format_cols.items()]

        self.header = line

        # Continue to read the file, this time the variants
        for j, line in enumerate(vcf):
            variant = Variant().process_variant(line, caller=self.caller)
            if variant.alt == '*':
                print("Warning: Vcf {} line {} has variant with alt=*".format(self.caller, str(i+j+1)))
            cleaned_variant = Variant.select_info(variant, info_cols, format_cols)
            # The dictionary is query by chr\tpos\tref\talt
            self.variants.update({cleaned_variant.variant_key: cleaned_variant})

        return self

    def process_somatic_vcf(self, cols, nid, tid):
        """Process somatic vcf, with normal and tumor sample id provided
        """

        vcf = open(self.name, 'r')
        info_dict, format_dict = od(), od()

        # Read the meta-information lines from the vcf
        for i, line in enumerate(vcf):
            # Handle exceptions: skip the AF and DP in INFO
            if (line.startswith('##INFO=<ID=DP') or
               line.startswith('##INFO=<ID=AF')):
                pass
            # Select the INFO/FORMAT lines
            elif line.startswith('##FORMAT'):
                vcf_header = VcfHeader(line)
                format_dict.update({vcf_header.meta_id: vcf_header})            
            elif line.startswith('##INFO'):
                vcf_header = VcfHeader(line)
                info_dict.update({vcf_header.meta_id: vcf_header})
            # Keep other meta-info lines
            elif line.startswith('##'):
                if line.startswith('##source='):
                    self.caller = line.replace('##source=', '').strip()
                self.meta_info.append(line)
            else:
                break

        if not self.caller:
            sys.exit("Cannot locate caller from file {}\nPlease add caller \
                     identify line '##source=(caller name)' to vcf header"
                     .format(self.name))

        # Select the columns from INFO/FORMAT
        info_cols, format_cols = extract_cols_somatic(info_dict, format_dict,
                                                      cols)

        # Add the INFO line (with caller) / FORMAT (unchanged) to header_list
        self.meta_info += [VcfHeader.write(VcfHeader.add_caller(v,
                           self.caller)) for k, v in info_cols.items()]
        self.meta_info += [VcfHeader.write(v) for k, v in format_cols.items()]

        self.header = line

        # find the location of normal / tumor column index
        if self.caller == 'strelka':
            normal_index, tumor_index = (self.header.split().index('NORMAL'),
                                         self.header.split().index('TUMOR'))
        elif (self.header.split()[9] == nid) and \
             (self.header.split()[10] == tid):
            normal_index, tumor_index = 9, 10
        elif (self.header.split()[10] == nid) and \
             (self.header.split()[9] == tid):
            normal_index, tumor_index = 10, 9
        else:
            sys.exit("Normal sample id [{}] or tumor sample id [{}] didn't match with file {}: [{}], [{}]"
                     .format(nid, tid, self.name, self.header.split()[9], self.header.split()[10]))

        # Continue to read the file, this time the variants
        for j, line in enumerate(vcf):
            variant = Variant().process_somatic_variant(
                    line, self.caller, normal_index, tumor_index)
            if variant.alt == '*':
                print("Warning: Line {} contains variant with alt=*".format(str(i+j+1)))
            cleaned_variant = Variant.select_info(variant, info_cols, format_cols,
                                                  caller=self.caller, somatic=True)
            # The dictionary is query by chr\tpos\tref\talt
            self.variants.update({cleaned_variant.variant_key: cleaned_variant})

        return self


