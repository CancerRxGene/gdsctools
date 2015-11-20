import pandas as pd
from gdsctools.report import HTMLTable


def test_htmltable():
    df = pd.DataFrame({
        'A':[1,2,10], 
        'B':[1,10,2], 
        'C':[1,10,2], 
        'url':['A', 'B', 'C']})
    html = HTMLTable(df, 'test')
    html.add_href('url')
    html.add_bgcolor('A')
    html.add_bgcolor('B', mode='clip', threshold=2)
    html.add_bgcolor('C', mode='max', threshold=2)
    print(html.to_html())


