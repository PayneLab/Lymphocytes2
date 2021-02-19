import pandas as pd
import requests
import os.path
import bs4
import requests
import urllib3
import csv
from os import path

#Data loader functions belong here. This is where
#  information about the data files is found.

def load_proteomics(version='current', level='protein', 
                    prefix="", suffix="Total Intensity",
                    contains=[], prepend_label=""):
    
    # Step 1. Which file are we reading from?
    file = get_file(key = version)
    if file==1:
        print("Error with file download.")
        return False
    
    # Step 2. Read the file    
    df = pd.read_csv(file, sep='\t', header=0, index_col=index_col)
    
    #We want to use "Protein" "Protein IDs" or "Protein ID" as the index, 
    #    or "Sequence" in the case of peptides;
    #    here we check which exists in this file.
    index = False
    if level == 'peptide': index = "Sequence"
    elif level == "protein":
        if df.columns.contains('Protein'): index = "Protein"
        elif df.columns.contains("Protein IDs"): index = "Protein IDs"
    if index:
        df.set_index(index)
    
    # Step 3. Filter 
    headings = df.columns
    
    if suffix:#filter by options such as suffix
        headings = [i for i in headings if i.endswith(suffix)]
    if prefix:#filter by columns beginning in prefix
        headings = [i for i in headings if i.startswith(prefix)]
    for req in contains:
        headings = [i for i in headings if req in i]
    for req in not_contains:
        headings = [i for i in headings if req not in i]
    
    # Optional 3b: Filter rows
    #    This may not be necessary or may be done differently.
    #    For example, ignoring those MQ marks as likely contaminants
    
    #drop contaminents and decoys
    if 'Potential contaminant' in df.headers:
        df = df.drop(df[df['Potential contaminant'] == '+'].index)
    if 'Reverse' in df.headers:
        df = df.drop(df[df.Reverse == '+'].index)
    if level=='protein':
        #optionally, discard those that were only identified by site
        #this will not work for peptide level analysis
        if 'Only identified by site' in df.headers:
            df = df.drop(df[df['Only identified by site'] == '+'].index)
    
    df = df[headings]
    
    # Step 4. Clean headers
    # Remove the prefix (ie, "Total Intensity") from the column names
    # optionally prepends a sample type (ie, "HeLa")
    new_names={}
    for c in df.columns.values: 
        sample_name = c[len(prefix):].strip()
        sample_name = c[:len(suffix)].strip()
        new_names[c] = "{0}_{1}".format(prepend_label, sample_name)
    df.rename(columns=new_names, inplace=True)
    df.head()
    
    # Return data
    return df
    
    

def load_max_quant(version = 'current', level='protein', 
                   prefix="Intensity", contains=["_"],
                   sample_type=""
                  ):
    #Takes a file and returns a dataframe.
    #    file: the file path to read from
    #    The rest of the paramters are used to select the columns.
    #    By default, it will look for ones starting with 'Reporter intensity'
    #        that do not contain 'count' or 'corrected' and use the 'Protein IDs'
    #        column as the indecies. These will be the raw intensity values.
    #file = get_file(key = version)#We need to add max_quant files to the index_url on box so we can use their keys on this
        
    if level=='protein':
        path = "data/proteinGroups_{0}.txt".format(version)
        url = "data/proteinGroups_{0}_url.txt".format(version)
    elif level=="peptide":
        path = "data/peptides_{0}.txt".format(version)
        url = "data/peptides_{0}_url.txt".format(version)
    else:
        #unknown level
        print ("Please specify either 'protein' or 'peptide' level.")
        return False
    file = download_file(download_to_path=path, url_file_path=url)

    
    #read in data
    df = pd.read_csv(file, sep='\t', header=0, index_col=0)
    
    #filter the columns based on the prefix and other "contains" requirements
    headings = df.columns
    if prefix:#filter by columns beginning in prefix
        headings = [i for i in headings if i.startswith(prefix)]
    for req in contains:
        headings = [i for i in headings if req in i]
        
    #drop contaminents and decoys
    df = df.drop(df[df['Potential contaminant'] == '+'].index)
    df = df.drop(df[df.Reverse == '+'].index)
    
    if level=='protein':
        #optionally, discard those that were only identified by site
        #this will not work for peptide
        df = df.drop(df[df['Only identified by site'] == '+'].index)
    
    df = df[headings]
    
    # Remove the prefix (ie, "Total Intensity") from the column names
    # optionally prepends a sample type (ie, "HeLa"
    new_names={}
    for c in df.columns.values: 
        sample_name = c[len(prefix):].strip()
        new_names[c] = "{0}_{1}".format(sample_type, sample_name)
    df.rename(columns=new_names, inplace=True)
    df.head()

    return df

def get_file(key = 'current'):
    #Takes the version we are looking for and sets up a table
    #from the url file so that we can use the version passed in as
    #a key to identify what url from the index table to download.
    url_file = open('data/index_url.txt', 'r')
    url = url_file.read().strip()
    url_file.close()
    
    table_file_path = download_file(download_to_path="data/index_table.tsv", url = url)
    table = pd.read_csv(table_file_path, sep='\t', header = 0, index_col = 'key')
    file_url = table.loc[key]
    
    file_name="data/{0}.tsv".format(key)
    url_name = file_url[0]
    
    
    return download_file(download_to_path=file_name, url=url_name, redownload = False)

def load_FragPipe(version = 'current', contains=[],level='protein', 
    suffix="Total Intensity"):
    #Takes a file and returns a dataframe.
    #    file: the file path to read from
    #    The rest of the paramters are used to select the columns.
    #    By default, it will look for ones ending with 'Total intensity'
    #        that do not contain 'count' or 'corrected' and use the 'Protein IDs'
    #        column as the indecies. These will be the raw intensity values.
    file = get_file(key = version)
    if file==1:
        print("Error with file download.")
        return False
        
    if version=='June':not_contains=['15']#drop extra replicate - Yiran said these two weren't good quality, I just forgot to not run it so for now I'll exclude it at this level
    else: not_contains=[]

        #read in data
    if level == 'protein': index_col = 3
    else: index_col=0 #for peptides and by default, take the first column as index
    df = pd.read_csv(file, sep='\t', header=0, index_col=index_col)
    
    #filter the columns based on the prefix and other "contains" requirements
    headings = df.columns
    
    if suffix:#filter by options such as suffix, contains
        headings = [i for i in headings if i.endswith(suffix)]
    for req in contains:
        headings = [i for i in headings if req in i]
    for req in not_contains:
        headings = [i for i in headings if req not in i]
    
    df = df[headings]
    
    # Remove the "Total Intensity" part of the column names
    new_names={}
    for c in df.columns.values: 
        new_names[c] = c.split(' ')[0]
    df.rename(columns=new_names, inplace=True)
    df.head()

    return df

def download_file(download_to_path="data/datafile.txt", url='', 
                  password_file_path="data/password.txt", redownload=False):
    """Download a file from a given url to the specified location.
    Parameters:
    path (str): The path to the file to save the file to on the local machine.
    Returns:
    str: The path the file was downloaded to.
    """
        
    if redownload or path.exists(download_to_path) == False: #If the file has been downloaded, or the user wants to update, download the file
        if url == '':
            print("URL MUST BE SPECIFIED FOR DOWNLOAD")
            return 1
        
        for i in range(2):

            with requests.Session() as session: # Use a session object to save cookies
                # Construct the urls for our GET and POST requests
                get_url = url
                post_url = get_url.replace("https://byu.box.com/shared", "https://byu.app.box.com/public")

                # Send initial GET request and parse the request token out of the response
                get_response = session.get(get_url)
                soup = bs4.BeautifulSoup(get_response.text, "html.parser")
                token_tag = soup.find(id="request_token")
                #print (token_tag)
                #print (type(token_tag))
                
                #This cheks if there is a password file and if it found a password requirement on the file
                if token_tag is not None:
                    #This identifies if the error was with the password file path.
                    if path.exists(password_file_path) == False:
                        print("MISSING PASSWORD FILE")
                        return 1
                
                    #print("Checking password...")
                    password_file = open(password_file_path, 'r')
                    password = password_file.read().strip()
                    password_file.close()
                    token = token_tag.get("value")

                    # Send a POST request, with the password and token, to get the data
                    payload = {
                        'password': password,
                        'request_token': token}
                    response = session.post(post_url, data=payload)

                    with open(download_to_path, 'wb') as dest:
                        dest.write(response.content)
                        
                
                #This will download the file if it was not password protected
                else:
                    #print("No password needed")
                    response = requests.get(post_url, allow_redirects=True)
                    with open(download_to_path, 'wb') as out_file:
                        out_file.write(response.content)
                
                    
    return download_to_path

def load_fasta(file="data/uniprot-filtered-proteome_3AUP000005640_reviewed_human.fasta"):
    
    #file is formated:
    #>sp|Q96IY4|CBPB2_HUMAN Carboxypeptidase B2 OS=Homo sapiens OX=9606 GN=CPB2 PE=1 SV=2
    #MKLCS...
    headings = {}
    with open(file) as f:
        for line in f:
            if line.startswith('>'):#header line
                ID = line.split('|')[1]
                name=line.split('|')[2].split('=')[0].strip('OS')
                headings[ID]=name
    headings = pd.Series(list(headings.values()), index=headings.keys())
    
    return headings

        
def names_max_quant():
    file = download_file(download_to_path="data/proteinGroups.txt", url_file_path="data/proteinGroups_url.txt")
    df = pd.read_csv(file, sep='\t', header=0, index_col=0, usecols=['Protein IDs','Gene names','Fasta headers'])

    return df
    
    
def names_FragPipe(month='June', contains=['Subject1']):
    file_name="data/combined_protein_{0}_FP.tsv".format(month)
    url_file_path="data/combined_protein_{0}_FP_url.txt".format(month)
    file = download_file(download_to_path=file_name, url_file_path=url_file_path)
    df = pd.read_csv(file, sep='\t', header=0, index_col=0, usecols=['Protein ID','Gene Names','Description'])
    return df             




#the following are for CBC and metabolics data
def to_num(x, exclude=['date', 'Testing entity']):
    if x.name not in exclude:
        x = pd.to_numeric(x)
    return x

def load_cbc(sub=0, time = '2020_09'):
    cbc_data = pd.concat([pd.read_csv("data/{0}/cbc_sub1.tsv".format(time), sep='\t'),
                          pd.read_csv("data/{0}/cbc_sub2.tsv".format(time), sep='\t')])
    units = cbc_data.iloc[0]
    cbc_data = cbc_data.drop(cbc_data['subject']==np.nan)
    cbc_data = cbc_data.apply(lambda x: to_num(x))
 
    cbc_data['date']=pd.to_datetime(cbc_data['date'])
    
    if sub: # if the user specifies they only want a particular subject
        cbc_data=cbc_data[cbc_data.apply(lambda x: x['subject']==sub, axis='columns')]
    return cbc_data


def load_metabolites(sub=0, time = '2020_09'):
    metab_data = pd.concat([pd.read_csv("data/{0}/metab_sub1.tsv".format(time), sep='\t'),
                            pd.read_csv("data/{0}/metab_sub2.tsv".format(time), sep='\t')], sort=False)
    metab_data['date']=pd.to_datetime(metab_data['date'])
    if sub:
        metab_data=metab_data[metab_data.apply(lambda x: x['Subject']==sub, axis='columns')]
    return metab_data

def load_reference(time = '2020_09'):
    met_ref =pd.read_csv("data/{0}/metab_reference.tsv".format(time), sep='\t', index_col=0)
    cbc_ref =pd.read_csv("data/{0}/cbc_reference.tsv".format(time), sep='\t', index_col=0)
    ref = pd.concat([met_ref, cbc_ref], axis=1)
    return ref
