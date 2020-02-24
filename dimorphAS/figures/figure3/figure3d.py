
import pandas as pd


from plotly.offline import iplot, init_notebook_mode
import plotly.offline as offline
import plotly.graph_objs as go
import plotly.io as pio

import os
import numpy as np

#init_notebook_mode(connected=True)


data3d = pd.read_csv('3dtable.txt',sep='\t',dtype='object')


Male = go.Scatter3d(
    x = data3d.loc[data3d['Sex']=='male',data3d.columns[0]],
    y = data3d.loc[data3d['Sex']=='male',data3d.columns[1]],
    z=  data3d.loc[data3d['Sex']=='male',data3d.columns[2]],
    mode = 'markers',name='Male', marker=dict(
        size=12,color='rgba(60,84,136,1)',

        opacity=0.8
    )
)

Female = go.Scatter3d(
    x = data3d.loc[data3d['Sex']=='female',data3d.columns[0]],
    y = data3d.loc[data3d['Sex']=='female',data3d.columns[1]],
    z=  data3d.loc[data3d['Sex']=='female',data3d.columns[2]],
    mode = 'markers',name='Female', marker=dict(
        size=12,color='rgba(220,0,0,1)',
        opacity=0.8
    )
)


layout = go.Layout(
    autosize=True,
    font=dict(family='Courier New, monospace',size=18,color='black'),
    scene = dict(
        xaxis = dict(
            title=list(data3d.columns.values)[0]),
        yaxis = dict(
            title=list(data3d.columns.values)[1]),
        zaxis = dict(
            title=list(data3d.columns.values)[2]),),
    width=2000,
    height=1000,
    legend=dict(x=0.45,y=0.3,font=dict(
            family='sans-serif',
            size=36)
    ),
    margin=dict(r=50, b=25,l=25, t=30),
)

data = [Male,Female]

fig = go.Figure(data=data, layout=layout)

offline.plot(fig, filename='basic-scatter')

#init_notebook_mode(connected=True)
#pio.write_image(fig, 'fig3d.pdf')

pio.write_image(fig, '3d-scatter-colorscale.pdf')
