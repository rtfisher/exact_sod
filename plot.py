import numpy as np
import plotly.express as px   #plotly is interactive plotting library. Supports hover function natively.
from plotly.subplots import make_subplots
import plotly.graph_objects as go

def one_plotter(filename1):
    fname1 = filename1 
    
    target = open(fname1,'r')
    lists = []
    var1 = []
    var2 = []
    var3 = []
    var4 = []
    

    i=0
    for line in target:
        arr = line.strip()
        data = arr.split()
        #print(data)
        lists.append([])
        for j in range(len(data)):
           lists[i].append(data[j])
        i=i+1

    #print list from exact_sod.out
    count = 0
    for e in lists:
        if isinstance(e, list):
            count = count + 1
    
    count = count - 1

    
    for i in range(count): 
        var1.append(float(lists[i][0])+0.5)
    
    for i in range(count): 
        var2.append(float(lists[i][1]))
        
    for i in range(count):  
        var3.append(float(lists[i][3]))
        
    for i in range(count):
        var4.append(float(lists[i][2]))
        
    
    dens = np.array(var2) # making the list to a proper numpy array so that it drops dimensions if any.
    pres = np.array(var3)
    velx = np.array(var4)
    xpos = np.array(var1)

    fig = make_subplots(rows=3, cols=1,shared_xaxes = False,
                    vertical_spacing=0.08)

    fig.append_trace(go.Scatter(
        x=xpos,
        y=dens,
        name = "Density"
    ), row=1, col=1)  
    
    fig.append_trace(go.Scatter(
        x=xpos,
        y=pres,
        name = "Pressure"
    ), row=2, col=1) 
    
    fig.append_trace(go.Scatter(
        x=xpos,
        y=velx,
        name = "x velocity"
    ), row=3, col=1)
    
    
    # Update xaxis properties
    fig.update_xaxes(title_text="X", row=1, col=1)
    fig.update_xaxes(title_text="X", row=2, col=1)
    fig.update_xaxes(title_text="X", row=3, col=1)
    
    # Update yaxis properties
    fig.update_yaxes(title_text="Density", row=1, col=1)
    fig.update_yaxes(title_text="Pressure", row=2, col=1)
    fig.update_yaxes(title_text="Vel X", row=3, col=1)
   
    fig.update_layout(height=800, width=800)
    fig.show()  # change to save if you want to


one_plotter("exact_sod.out") #make sure to include the correct path to this output file.