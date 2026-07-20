#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', 'Setting_Env_Variables_p2.py')
google_project_id = get_ipython().run_line_magic('env', 'GOOGLE_CLOUD_PROJECT')
bucket = os.getenv("WORKSPACE_BUCKET")

import os
import subprocess
import pandas as pd
from datetime import datetime

import hail as hl
hl.init(gcs_requester_pays_configuration=google_project_id) #, log=f'{bucket}/hail_logs')
hl.default_reference(new_default_reference = "GRCh38")


# In[2]:


auxiliary_path = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux"
print(f'aux path: {auxiliary_path}')
vat_path = f'{auxiliary_path}/vat/*.gz'
print(f'vat path: {vat_path}')

# mt_wgs_path = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt"
# flagged_samples = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv"
# array_path = "gs://vwb-aou-datasets-controlled/v9/microarray/hail.mt"
# # this path is in the official ct v8 fw
# bed_path = "gs://aou-tutorial-notebooks-wb-sunny-radish-6214/genomic_test_data/random_sites_GRCh38.txt"


# In[4]:


get_ipython().run_cell_magic('time', '', 'vat_table = hl.import_table(vat_path, force=True, quote=\'"\', delimiter="\\t", force_bgz=True,) # impute=True)\nvat_table = vat_table.add_index(name=\'id\')\n')


# In[6]:


gene = 'LMNA'
ens_id = 'ENSG00000160789'

filt_vat = vat_table.filter(vat_table["gene_id"]==ens_id)


# In[7]:


filt_vat.write(f'{bucket}/hail_files/lmna_vat.ht', overwrite=True)


# In[ ]:




