# analyze all models vs common10

import os
import pandas as pd
from scipy import stats

c10dir = '/where/common10/files/are/'

# initiate a dataframe to store results
comb_results = pd.DataFrame(columns=['test', 'model', 'rep','r95', 'psigcount','adjpsigcount'])

# Iterate through numbers 1 to 10
for i in range(1, 11):
    # Construct the filename
    filename = f'v43_100k_{i}.csv'
    # Full path to the file
    filepath = os.path.join(c10dir, filename)
    # Check if the file exists
    if os.path.exists(filepath):
        # Read the file into a dataframe
        df = pd.read_csv(filepath)
        
        distname = df.iloc[0,0]
        params = df.iloc[0,2]
        params = eval(params)
        #Load distribution
        dist = getattr(stats, distname)
        distribution = dist(*params)

        # get the 95% percentile which is the p=0.05 r value
        r95 = distribution.ppf(0.95)
    
    pfile=f'v43vXIST_100k_{i}.csv'
    pfilepath = os.path.join(c10dir, pfile)
    if os.path.exists(pfilepath):
        pdf = pd.read_csv(pfilepath)
        # rename the columns
        pdf.columns = ['gene', 'p']

        # get the p values
        p = pdf['p']
        # get the number of significant r values
        psigcount = sum(p<0.05)

    adjpfile=f'v43vXIST_100k_bh_{i}.csv'
    adjpfilepath = os.path.join(c10dir, adjpfile)
    if os.path.exists(adjpfilepath):
        adjpdf = pd.read_csv(adjpfilepath)
        # rename the columns
        adjpdf.columns = ['gene', 'adjp']

        # get the adjp values
        adjp = adjpdf['adjp']
        # get the number of significant r values
        adjpsigcount = sum(adjp<0.05)

    # add the results to the dataframe
    new_row={'test':'common10','model': distname, 'rep': i, 'r95': r95, 'psigcount': psigcount, 'adjpsigcount': adjpsigcount}
    new_row_df = pd.DataFrame([new_row])
    comb_results = pd.concat([comb_results, new_row_df], ignore_index=True)



# analyze all data

alldir = '/where/all/files/are/'

# Iterate through numbers 1 to 10
for i in range(1, 11):
    # Construct the filename
    filename = f'v43_100k_{i}.csv'
    # Full path to the file
    filepath = os.path.join(alldir, filename)
    # Check if the file exists
    if os.path.exists(filepath):
        # Read the file into a dataframe
        df = pd.read_csv(filepath)
        
        distname = df.iloc[0,0]
        params = df.iloc[0,2]
        params = eval(params)
        #Load distribution
        dist = getattr(stats, distname)
        distribution = dist(*params)

        # get the 95% percentile which is the p=0.05 r value
        r95 = distribution.ppf(0.95)
    
    pfile=f'v43vXIST_100k_{i}.csv'
    pfilepath = os.path.join(alldir, pfile)
    if os.path.exists(pfilepath):
        pdf = pd.read_csv(pfilepath)
        # rename the columns
        pdf.columns = ['gene', 'p']

        # get the p values
        p = pdf['p']
        # get the number of significant r values
        psigcount = sum(p<0.05)

    adjpfile=f'v43vXIST_100k_bh_{i}.csv'
    adjpfilepath = os.path.join(alldir, adjpfile)
    if os.path.exists(adjpfilepath):
        adjpdf = pd.read_csv(adjpfilepath)
        # rename the columns
        adjpdf.columns = ['gene', 'adjp']

        # get the adj p values
        adjp = adjpdf['adjp']
        # get the number of significant r values
        adjpsigcount = sum(adjp<0.05)

    # add the results to the dataframe
    new_row={'test':'all_bf1','model': distname, 'rep': i, 'r95': r95, 'psigcount': psigcount, 'adjpsigcount': adjpsigcount}
    new_row_df = pd.DataFrame([new_row])
    comb_results = pd.concat([comb_results, new_row_df], ignore_index=True)

    #############################
    # second best fit
    distname = df.iloc[1,0]
    params = df.iloc[1,2]
    params = eval(params)
    #Load distribution
    dist = getattr(stats, distname)
    distribution = dist(*params)
    # get the 95% percentile which is the p=0.05 r value
    r95 = distribution.ppf(0.95)
    
    pfile=f'v43vXIST_2_100k_{i}.csv'
    pfilepath = os.path.join(alldir, pfile)
    if os.path.exists(pfilepath):
        pdf = pd.read_csv(pfilepath)
        # rename the columns
        pdf.columns = ['gene', 'p']

        # get the p values
        p = pdf['p']
        # get the number of significant r values
        psigcount = sum(p<0.05)

    adjpfile=f'v43vXIST_2_100k_bh_{i}.csv'
    adjpfilepath = os.path.join(alldir, adjpfile)
    if os.path.exists(adjpfilepath):
        adjpdf = pd.read_csv(adjpfilepath)
        # rename the columns
        adjpdf.columns = ['gene', 'adjp']

        # get the adj p values
        adjp = adjpdf['adjp']
        # get the number of significant r values
        adjpsigcount = sum(adjp<0.05)

    # add the results to the dataframe
    new_row={'test':'all_bf2','model': distname, 'rep': i, 'r95': r95, 'psigcount': psigcount, 'adjpsigcount': adjpsigcount}
    new_row_df = pd.DataFrame([new_row])
    comb_results = pd.concat([comb_results, new_row_df], ignore_index=True)


    #############################
    # third best fit
    distname = df.iloc[2,0]
    params = df.iloc[2,2]
    params = eval(params)
    #Load distribution
    dist = getattr(stats, distname)
    distribution = dist(*params)
    # get the 95% percentile which is the p=0.05 r value
    r95 = distribution.ppf(0.95)
    
    pfile=f'v43vXIST_3_100k_{i}.csv'
    pfilepath = os.path.join(alldir, pfile)
    if os.path.exists(pfilepath):
        pdf = pd.read_csv(pfilepath)
        # rename the columns
        pdf.columns = ['gene', 'p']

        # get the p values
        p = pdf['p']
        # get the number of significant r values
        psigcount = sum(p<0.05)

    adjpfile=f'v43vXIST_3_100k_bh_{i}.csv'
    adjpfilepath = os.path.join(alldir, adjpfile)
    if os.path.exists(adjpfilepath):
        adjpdf = pd.read_csv(adjpfilepath)
        # rename the columns
        adjpdf.columns = ['gene', 'adjp']

        # get the adj p values
        adjp = adjpdf['adjp']
        # get the number of significant r values
        adjpsigcount = sum(adjp<0.05)

    # add the results to the dataframe
    new_row={'test':'all_bf3','model': distname, 'rep': i, 'r95': r95, 'psigcount': psigcount, 'adjpsigcount': adjpsigcount}
    new_row_df = pd.DataFrame([new_row])
    comb_results = pd.concat([comb_results, new_row_df], ignore_index=True)


# save the results to a csv file
comb_results.to_csv('allvscommon10.csv', index=False)
        
