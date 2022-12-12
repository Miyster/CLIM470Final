# Use matplotlib and cartopy to make plots for wind vectors, see paper for plot ideas (taylor will look into this once code is written and able to be exported)
# taylor idea/research dump:
  # https://scitools.org.uk/iris/docs/v2.2/examples/Meteorology/wind_speed.html
  # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.quiver.html
  # https://researchworkspace.com/file/2588789/plot-wind-vectors-from-grid.ipynb?preview=true
  # find way to create wind barbs

  
  
#clim 470 final project

# d = 5e+05

  #delt = 5*60
  #d = 5e+05

    import matplotlib.pyplot as plt
    import numpy as np
    import os
    import pandas as pd
    import numpy as np
    import time
    import csv

    path='/homes/tvineyar/miocenedata/'
    delt_file=pd.read_csv(path+"delt_5_60_d_5e_05.csv",header=0,error_bad_lines=False, sep=',')
    delt_file

    data0=delt_file.values[:,:]
    data0

    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    x,y = np.meshgrid(np.linspace(-10,10,25),np.linspace(-10,10,25))

    u = np.asarray(data0[:,1],dtype=np.float64)
    v = np.asarray(data0[:,2],dtype=np.float64)

    plt.title('delta t = 5*60, d =5e05', fontsize=10)
    plt.quiver(x,y,u,v)
    plt.show()

  # delta t 10*60
  # d = 5e05

    path='/homes/tvineyar/miocenedata/'
    delt_file=pd.read_csv(path+"deltat_10_60_d_5e_05.csv",header=0,error_bad_lines=False, sep=',')
    delt_file

    data0=delt_file.values[:,:]
    data0

    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    x,y = np.meshgrid(np.linspace(-10,10,15),np.linspace(-10,10,15))

    u = np.asarray(data0[:,1],dtype=np.float64)
    v = np.asarray(data0[:,2],dtype=np.float64)

    plt.title('delta t = 10*60, d =5e05', fontsize=10)
    plt.quiver(x,y,u,v)
    plt.show()

  # delta t = 2.5*60
  # d = 5e+05

    path='/homes/tvineyar/miocenedata/'
    delt_file=pd.read_csv(path+"delt_2.5_60_d_5e_05.csv",header=0,error_bad_lines=False, sep=',')
    delt_file

    data0=delt_file.values[:,:]
    data0

    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    x,y = np.meshgrid(np.linspace(-10,10,35),np.linspace(-10,10,35))

    plt.title('delta t = 2.5*60, d =5e05', fontsize=10)
    u = np.asarray(data0[:,1],dtype=np.float64)
    v = np.asarray(data0[:,2],dtype=np.float64)

    plt.quiver(x,y,u,v)
    plt.show()

# d = 2.5e+05


  # delta t = 10*60
  # d = 2.5e+0.5
  
    path='/homes/tvineyar/miocenedata/'
    delt_file=pd.read_csv(path+"delt_10_60_d_2.5_0.5.csv",header=0,error_bad_lines=False, sep=',')
    delt_file
    
    data0=delt_file.values[:,:]
    data0
    
    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    x,y = np.meshgrid(np.linspace(-10,10,15),np.linspace(-10,10,15))

    plt.title('delta t = 10*60, d =2.5e+05', fontsize=10)
    u = np.asarray(data0[:,1],dtype=np.float64)
    v = np.asarray(data0[:,2],dtype=np.float64)

    plt.quiver(x,y,u,v)
    plt.show()
    
  # delta t = 5*60
  # d = 2.5e+05
  
    path='/homes/tvineyar/miocenedata/'
    delt_file=pd.read_csv(path+"delt_5_60_d_2.5e_05.csv",header=0,error_bad_lines=False, sep=',')
    delt_file 
    
    data0=delt_file.values[:,:]
    data0
    
    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    x,y = np.meshgrid(np.linspace(-10,10,25),np.linspace(-10,10,25))

    plt.title('delta t = 5*60, d =2.5e+05', fontsize=10)
    u = np.asarray(data0[:,1],dtype=np.float64)
    v = np.asarray(data0[:,2],dtype=np.float64)

    plt.quiver(x,y,u,v)
    plt.show()

  # delta t = 2.5*10
  # d = 2.5e+0.5
  
    path='/homes/tvineyar/miocenedata/'
    delt_file=pd.read_csv(path+"delt_2.5_60_d_2.5e_05.csv",header=0,error_bad_lines=False, sep=',')
    delt_file
    
    data0=delt_file.values[:,:]
    data0
    
    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    x,y = np.meshgrid(np.linspace(-10,10,35),np.linspace(-10,10,35))

    plt.title('delta t = 2.5*60, d =2.5e+05', fontsize=10)
    u = np.asarray(data0[:,1],dtype=np.float64)
    v = np.asarray(data0[:,2],dtype=np.float64)

    plt.quiver(x,y,u,v)
    plt.show()
    
# d = 1.25e+05

  # delta t = 10*60
  # d = 1.25e+05
  
    path='/homes/tvineyar/miocenedata/'
    delt_file=pd.read_csv(path+"delt_10_60_d_1.25e_05.csv",header=0,error_bad_lines=False, sep=',')
    delt_file
    
    data0=delt_file.values[:,:]
    data0
    
    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    x,y = np.meshgrid(np.linspace(-10,10,15),np.linspace(-10,10,15))

    plt.title('delta t = 10*60, d =1.25e+05', fontsize=10)
    u = np.asarray(data0[:,1],dtype=np.float64)
    v = np.asarray(data0[:,2],dtype=np.float64)

    plt.quiver(x,y,u,v)
    plt.show()
    
  #delta t = 5*60
  #d = 1.25e+05
  
    path='/homes/tvineyar/miocenedata/'
    delt_file=pd.read_csv(path+"delt_5_60_d_1.25_05.csv",header=0,error_bad_lines=False, sep=',')
    delt_file
    
    data0=delt_file.values[:,:]
    data0

    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    x,y = np.meshgrid(np.linspace(-10,10,25),np.linspace(-10,10,25))

    plt.title('delta t = 5*60, d =1.25e+05', fontsize=10)
    u = np.asarray(data0[:,1],dtype=np.float64)
    v = np.asarray(data0[:,2],dtype=np.float64)

    plt.quiver(x,y,u,v)
    plt.show()
    
  # delta t = 2.5*60
  # d = 1.25e+05
  
    path='/homes/tvineyar/miocenedata/'
    delt_file=pd.read_csv(path+"delt_2.5_60_d_1.25e_05.csv",header=0,error_bad_lines=False, sep=',')
    delt_file
    
    data0=delt_file.values[:,:]
    data0
    
    import numpy as np
    import matplotlib.pyplot as plt
    %matplotlib inline

    x,y = np.meshgrid(np.linspace(-10,10,35),np.linspace(-10,10,35))

    plt.title('delta t = 2.5*60, d =1.25e+05', fontsize=10)
    u = np.asarray(data0[:,1],dtype=np.float64)
    v = np.asarray(data0[:,2],dtype=np.float64)

    plt.quiver(x,y,u,v)
    plt.show()
    
    

  


