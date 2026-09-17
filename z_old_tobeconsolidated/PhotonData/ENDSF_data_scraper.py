import pandas as pd
import urllib.request
import mendeleev as mdl
import time
from mendeleev.fetch import fetch_table

# the service URL
livechart = "https://nds.iaea.org/relnsd/v1/data?fields=gammas&nuclides="

isotopes = fetch_table("isotopes", index_col="id")


def lc_pd_dataframe(url, save_as=None):
    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
    response = urllib.request.urlopen(req)
    raw_data = response.read().decode('utf-8').strip()
    if raw_data != "0":
        with open(save_as, "w") as f:
            f.write(raw_data)
        print(f'written data for {url}')
    return

# i=10
# element_str = str(mdl.element(int(isotopes.at[i,'atomic_number'])).symbol).lower()
# mass_str = str(isotopes.at[i,'mass_number'])
# req_add = str(mass_str + element_str)
# req_add = '12c'
# print(livechart + req_add)
# lc_pd_dataframe(livechart + req_add + '&download=csv', save_as=f'.\\{req_add}.csv')

if __name__ == "__main__":
    for i in range(1,len(isotopes)):
        element_str = str(mdl.element(int(isotopes.at[i,'atomic_number'])).symbol).lower()
        mass_str = str(isotopes.at[i,'mass_number'])
        req_add = str(mass_str + element_str)
        lc_pd_dataframe(livechart + req_add, save_as=f'.\\gamma_data\\{req_add}.csv')
        time.sleep(0.5)