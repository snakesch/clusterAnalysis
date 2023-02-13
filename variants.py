#!/usr/bin/env python3

class Variant:
        
    def __init__(self, bp, p, chrom, snp_id):
        
        import numpy as np 
    
        self.bp = bp
        self.p = p
        self.chrom = chrom
        self.logP = -np.log10(p)
        self.snpid = snp_id
        self.r2 = -1
    
    def __eq__(self, other):
        
        return self.bp == other.bp and self.chrom == other.chrom and self.p == other.p
    
    def __repr__(self):
        '''
        String representation of variant
        '''
        if self.r2 != -1:
            return "chr" + str(self.chrom) + ":" + str(self.bp) + " " \
                   "-logP = " + str(self.logP) + "; r2 = " + str(self.r2) + "\n"
        else:
            return "chr" + str(self.chrom) + ":" + str(self.bp) + " " \
                   "-logP = " + str(self.logP) + "\n"
    
class Cluster:

    def __init__(self, var_lst: list[Variant]):
        self.variants = var_lst
        self.update_stat()
        self.chrom = var_lst[0].chrom
        self.ld_ref = None
        
    def __len__(self):
        '''
        Returns total number of variants in the clsuter
        '''
        return len(self.variants)
    
    def __repr__(self):
        
        return "chr" + str(self.chrom) + ":" + str(self.start) + "-" + str(self.end)
    
    def __eq__(self, other):
        
        return self.variants == other.variants and self.chrom == other.chrom

    @classmethod
    def add_variant(self, *var: Variant):
        for _var in var:
            if _var not in self.variants:
                self.variants.extend([_var])
        self.update_stat()

    @classmethod
    def update_stat(self):
        self.maxp = max([var.p for var in self.variants])
        self.minp = min([var.p for var in self.variants])
        self.start = min([var.bp for var in self.variants])
        self.end = max([var.bp for var in self.variants])
    
    @classmethod
    def get_significant_vars(self, p_thresh):
        '''
        Get variants with p < p_thresh
        '''
        sig_vars = []
        for var in self.variants:
            if var.p < p_thresh:
                sig_vars.append(var)
        
        return sig_vars
    
    @classmethod
    def set_ld_reference(self, var: Variant):
        '''
        Set LD reference SNP
        '''
        self.ld_ref = var
    
    @classmethod
    def add_r2_info(self, ld_df, ld_ref: Variant):
        
        import logging
        logger = logging.getLogger("root")
        
        self.ld_ref = ld_ref
        logger.debug(f"Setting LD refence at {self.ld_ref.snpid}")
        ld_df = ld_df[(ld_df["CHR_A"] == self.chrom) & (ld_df["BP_A"] == self.ld_ref.bp)]
        r2_dict = {k: v for k, v in zip(ld_df["BP_B"], ld_df["R2"])}
        for other_var in self.variants:
            other_var.r2 = r2_dict.get(other_var.bp, -1)
    
    @classmethod
    def clear_r2(self):
        '''
        Reset r2 info for each member of the cluster
        '''
        self.ld_ref = None
        for var in self.variants:
            var.r2 = -1

    @classmethod
    def get_outliers(self):
        
        import numpy as np
        
        logP = sorted([var.logP for var in self.variants])
        qts = np.percentile(logP, [25, 75])
        iqr = qts[1] - qts[0]
        lower = qts[0] - 1.5 * iqr
        upper = qts[1] + 1.5 * iqr
        
        n_outliers = 0
        for val in logP:
            if val < lower:
                n_outliers = n_outliers + 1
            else:
                break
        for val in logP[::-1]:
            if val > upper:
                n_outliers = n_outliers + 1
            else:
                break
        
        return n_outliers