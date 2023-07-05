from bs4 import BeautifulSoup
import requests 
import csv
import os

# Reading in a HTML from a website
url = "https://people.mbi.ucla.edu/sawaya/amyloidatlas/"
result = requests.get(url)
doc = BeautifulSoup(result.text, "html.parser")
#print(doc.prettify())

# Writing the contents of `doc` to a text file
#with open("AA_HTML.txt", "w") as f:
#    f.write(str(doc))

# Accessing Table Data = <td><big><centre> data of interest </centre></big></td> in html terms
table_values = doc.find_all("td")

# initializing headers
headers = []

# Initializing lists
row_data = []
current_row = []


# The first 9 values in table_values are the column headers
# Of which 6 and 7 we do not want (they contain images)
# for loop extracts information from every row except for index 5 and 6 because of negative indexing
   
for i in range(0,len(table_values)):
    if i < 9 and i % 9 != 5 and i % 9 != 6:
        headers.append(table_values[i].get_text(strip=True))
    
    if i >= 9 and i % 9 != 5 and i % 9 != 6:
        x = table_values[i].get_text(strip=True)

        current_row.append(x)
        #current_row.append(table_values[i].get_text(strip=True))

    if i > 8 and i % 9 == 8:
        row_data.append(current_row)
        current_row = []    

# adding headers to the beginning of the list
row_data.insert(0,headers) 


# Creating Output Directory if it does not exist
if not os.path.exists("Output"):
    os.makedirs("Output")

# Saving list of lists to a text file
with open("Output/atlas_scraping.txt", "w", newline="") as file:
    writer = csv.writer(file, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
    writer.writerows(row_data)