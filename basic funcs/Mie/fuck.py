#encoding=utf-8

import requests,re,os,codecs

url='http://www.iap.unibe.ch/publications/download/2002-11/'

r = requests.get(url)

for line in r.text.split('\n'):
    r = re.search(r'<a href="(.+?\.m)">',line)
    if r:
        a=r.group(1)
        sb=requests.get(url+a)
        with codecs.open(a,'w','utf-8') as f:
            f.write(sb.text)
        

