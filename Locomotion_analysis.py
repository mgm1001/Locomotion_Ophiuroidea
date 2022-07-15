#!/usr/bin/env python
# coding: utf-8

# ## Locomotion analysis of five-armed ophiuroids

# Reference: Goharimanesh, M., Stöhr, S., Ghassemzadeh, F., Mirshamsi, O., Adriaens, D. Arm kinematics in Ophiuroidea (Echinodermata): an analysis on pattern consistency in arm use across species. 

# In below, the calculation of four kinematic variables (V1: Sinuoisity, V2: Arm angle, V3: Disc displacement, V4: Slip angle) between frames 1-5 for arms 1-5 is shown.

# ### Frame 1

# In[ ]:


#######  To choose the cooridnates of Arm 1 (CURVE1) ########
    
x1_F1=df_1['X-Coordinate '] [df_1['Curve Name']=='CURVE 1']
y1_F1=df_1['Y-Coordinate '] [df_1['Curve Name']=='CURVE 1']

######  Axes transformation ########
x121=np.array((x1_F1-float(x1_F1[:1])))
y121=np.array((y1_F1-float(y1_F1[:1])))

######  Axes rotation  ########
m=(y121[-1])/x121[-1]
angle1= np.degrees(np.arctan(m))

rules1 = [x121[-1]>0, y121[-1]>0, m>0]
rules2 = [x121[-1]>0, y121[-1]<0, m<0]
rules3 = [x121[-1]<0, y121[-1]<0, m>0]
rules4 = [x121[-1]<0, y121[-1]>0, m<0]


if  all(rules1):
    theta1= 360-angle1
elif  all(rules2): 
    theta1= abs(angle1)
elif  all(rules3): 
    theta1= 180-angle1
elif  all(rules4): 
    theta1= 180-angle1
print(theta1)

x_F1=x121*(math.cos(math.radians(theta1)))-y121*(math.sin(math.radians(theta1)))
y_F1=x121*(math.sin(math.radians(theta1)))+y121*(math.cos(math.radians(theta1)))



# ## Variable 1: Sinuosity for each arm (1-5) within each time frame (1-5)

# The calculation of sinuosity in arms 2-5 is written in script below Frames 1-5.

# In[ ]:


######   Equidistant landmarks   ########
f=interp1d(x_F1,y_F1)
distance = np.cumsum(np.sqrt( np.ediff1d(x_F1, to_begin=0)**2 + np.ediff1d(y_F1, to_begin=0)**2 ))
distance = distance/distance[-1]

fx, fy = interp1d( distance, x_F1 ), interp1d( distance, y_F1 )

alpha = np.linspace(0, 1, 15)
x_regular_F1, y_regular_F1 = fx(alpha), fy(alpha)


################################to calculate sinuosity of arm 1 in Frame 1 ##################################

#calculating the length of a curve
d=np.zeros(14)
for i in range(14):
    d[i]=np.sqrt((x_regular_F1[i+1]-x_regular_F1[i])**2+(y_regular_F1[i+1]-y_regular_F1[i])**2)
print(d)
length=sum(d)

#calculating the distance between first and the last point (shortest length)
distance=np.sqrt((x_regular_F1[14]-x_regular_F1[0])**2+(y_regular_F1[14]-y_regular_F1[0])**2)

#calculating siniosity
sinuosity1_F1=length/distance


# In[ ]:



                                #######  To choose the cooridnates of Arm 2 (CURVE2) ########
    
x11_F1=df_1['X-Coordinate '] [df_1['Curve Name']=='CURVE 2']
y11_F1=df_1['Y-Coordinate '] [df_1['Curve Name']=='CURVE 2']
x122=np.array((x11_F1-float(x11_F1[:1])))
y122=np.array((y11_F1-float(y11_F1[:1])))

m2=(y122[-1])/x122[-1]
angle2 = np.degrees(np.arctan(m2))

rules1 = [x122[-1]>0, y122[-1]>0, m2>0]
rules2 = [x122[-1]>0, y122[-1]<0, m2<0]
rules3 = [x122[-1]<0, y122[-1]<0, m2>0]
rules4 = [x122[-1]<0, y122[-1]>0, m2<0]

if  all(rules1):
    theta2= 360-angle2
elif  all(rules2): 
    theta2= abs(angle2)
elif  all(rules3): 
    theta2= 180-angle2
elif  all(rules4): 
    theta2= 180-angle2
print(theta2)

x2_F1=x122*(math.cos(math.radians(theta2)))-y122*(math.sin(math.radians(theta2)))
y2_F1=x122*(math.sin(math.radians(theta2)))+y122*(math.cos(math.radians(theta2)))

f=interp1d(x2_F1,y2_F1)

distance = np.cumsum(np.sqrt( np.ediff1d(x2_F1, to_begin=0)**2 + np.ediff1d(y2_F1, to_begin=0)**2 ))
distance = distance/distance[-1]

fx2, fy2 = interp1d( distance, x2_F1 ), interp1d( distance, y2_F1 )

alpha = np.linspace(0, 1, 15)
x_regular2_F1, y_regular2_F1 = fx2(alpha), fy2(alpha)

#to calculate sinuosity
d2=np.zeros(14)
for i in range(14):
    d2[i]=np.sqrt((x_regular2_F1[i+1]-x_regular2_F1[i])**2+(y_regular2_F1[i+1]-y_regular2_F1[i])**2)
print(d2)

length2=sum(d2)
distance2=np.sqrt((x_regular2_F1[14]-x_regular2_F1[0])**2+(y_regular2_F1[14]-y_regular2_F1[0])**2)
sinuosity2_F1=length2/distance2


                             #######  To choose the cooridnates of Arm 3 (CURVE3) ########

x12_F1=df_1['X-Coordinate '] [df_1['Curve Name']=='CURVE 3']
y12_F1=df_1['Y-Coordinate '] [df_1['Curve Name']=='CURVE 3']
x123=np.array((x12_F1-float(x12_F1[:1])))
y123=np.array((y12_F1-float(y12_F1[:1])))

m3=(y123[-1])/x123[-1]
angle3 = np.degrees(np.arctan(m3))

rules1 = [x123[-1]>0, y123[-1]>0, m3>0]
rules2 = [x123[-1]>0, y123[-1]<0, m3<0]
rules3 = [x123[-1]<0, y123[-1]<0, m3>0]
rules4 = [x123[-1]<0, y123[-1]>0, m3<0]

if  all(rules1):
    theta3= 360-angle3
elif  all(rules2): 
    theta3= abs(angle3)
elif  all(rules3): 
    theta3= 180-angle3
elif  all(rules4): 
    theta3= 180-angle3
print(theta3)


x3_F1=x123*(math.cos(math.radians(theta3)))-y123*(math.sin(math.radians(theta3)))
y3_F1=x123*(math.sin(math.radians(theta3)))+y123*(math.cos(math.radians(theta3)))


f=interp1d(x3_F1,y3_F1)
distance = np.cumsum(np.sqrt( np.ediff1d(x3_F1, to_begin=0)**2 + np.ediff1d(y3_F1, to_begin=0)**2 ))
distance = distance/distance[-1]

fx3, fy3 = interp1d( distance, x3_F1 ), interp1d( distance, y3_F1 )

alpha = np.linspace(0, 1, 15)
x_regular3_F1, y_regular3_F1 = fx3(alpha), fy3(alpha)


#to calculate sinuosity
d3=np.zeros(14)
for i in range(14):
    d3[i]=np.sqrt((x_regular3_F1[i+1]-x_regular3_F1[i])**2+(y_regular3_F1[i+1]-y_regular3_F1[i])**2)
print(d3)

length3=sum(d3)
distance3=np.sqrt((x_regular3_F1[14]-x_regular3_F1[0])**2+(y_regular3_F1[14]-y_regular3_F1[0])**2)
sinuosity3_F1=length3/distance3


                           #######  To choose the cooridnates of Arm 4 (CURVE4) ########
    
x13_F1=df_1['X-Coordinate '] [df_1['Curve Name']=='CURVE 4']
y13_F1=df_1['Y-Coordinate '] [df_1['Curve Name']=='CURVE 4']
x124=np.array((x13_F1-float(x13_F1[:1])))
y124=np.array((y13_F1-float(y13_F1[:1])))

m4=(y124[-1])/x124[-1]
angle4 = np.degrees(np.arctan(m4))

rules1 = [x124[-1]>0, y124[-1]>0, m4>0]
rules2 = [x124[-1]>0, y124[-1]<0, m4<0]
rules3 = [x124[-1]<0, y124[-1]<0, m4>0]
rules4 = [x124[-1]<0, y124[-1]>0, m4<0]
if  all(rules1):
    theta4= 360-angle4
elif  all(rules2): 
    theta4= abs(angle4)
elif  all(rules3): 
    theta4= 180-angle4
elif  all(rules4): 
    theta4= 180-angle4
print(theta4)

x4_F1=x124*(math.cos(math.radians(theta4)))-y124*(math.sin(math.radians(theta4)))
y4_F1=x124*(math.sin(math.radians(theta4)))+y124*(math.cos(math.radians(theta4)))


f=interp1d(x4_F1,y4_F1)
distance = np.cumsum(np.sqrt( np.ediff1d(x4_F1, to_begin=0)**2 + np.ediff1d(y4_F1, to_begin=0)**2 ))
distance = distance/distance[-1]

fx4, fy4 = interp1d( distance, x4_F1 ), interp1d( distance, y4_F1 )

alpha = np.linspace(0, 1, 15)
x_regular4_F1, y_regular4_F1 = fx4(alpha), fy4(alpha)

#to calculate sinuosity
d4=np.zeros(15)
for i in range(14):
    d4[i]=np.sqrt((x_regular4_F1[i+1]-x_regular4_F1[i])**2+(y_regular4_F1[i+1]-y_regular4_F1[i])**2)
print(d4)

length4=sum(d4)
distance4=np.sqrt((x_regular4_F1[14]-x_regular4_F1[0])**2+(y_regular4_F1[14]-y_regular4_F1[0])**2)
sinuosity4_F1=length4/distance4


                         #######  To choose the cooridnates of Arm 5 (CURVE5) ########
    
x14_F1=df_1['X-Coordinate '] [df_1['Curve Name']=='CURVE 5']
y14_F1=df_1['Y-Coordinate '] [df_1['Curve Name']=='CURVE 5']
x125=np.array((x14_F1-float(x14_F1[:1])))
y125=np.array((y14_F1-float(y14_F1[:1])))

m5=(y125[-1])/x125[-1]
angle5= np.degrees(np.arctan(m5))

rules1 = [x125[-1]>0, y125[-1]>0, m5>0]
rules2 = [x125[-1]>0, y125[-1]<0, m5<0]
rules3 = [x125[-1]<0, y125[-1]<0, m5>0]
rules4 = [x125[-1]<0, y125[-1]>0, m5<0]

if  all(rules1):
    theta5= 360-angle5
elif  all(rules2): 
    theta5= abs(angle5)
elif  all(rules3): 
    theta5= 180-angle5
elif  all(rules4): 
    theta5= 180-angle5
print(theta5)

x5_F1=x125*(math.cos(math.radians(theta5)))-y125*(math.sin(math.radians(theta5)))
y5_F1=x125*(math.sin(math.radians(theta5)))+y125*(math.cos(math.radians(theta5)))

f=interp1d(x5_F1,y5_F1)
distance = np.cumsum(np.sqrt( np.ediff1d(x5_F1, to_begin=0)**2 + np.ediff1d(y5_F1, to_begin=0)**2 ))
distance = distance/distance[-1]

fx5, fy5 = interp1d( distance, x5_F1), interp1d( distance, y5_F1 )

alpha = np.linspace(0, 1, 15)
x_regular5_F1, y_regular5_F1 = fx5(alpha), fy5(alpha)


#to calculate sinuosity
d5=np.zeros(14)
for i in range(14):
    d5[i]=np.sqrt((x_regular5_F1[i+1]-x_regular5_F1[i])**2+(y_regular5_F1[i+1]-y_regular5_F1[i])**2)
print(d5)

length5=sum(d5)
distance5=np.sqrt((x_regular5_F1[14]-x_regular5_F1[0])**2+(y_regular5_F1[14]-y_regular5_F1[0])**2)
sinuosity5_F1=length5/distance5


# ## Variable 2: Arm angle in each frame

# The scipts of calculating arm angle in frames  2-5 are written unders the specific time frame codes.

# In[ ]:


########## Arm angle F1 ############

###### To calculate the disc center point in F1#######
x_all=np.zeros(5)
x_all[0]=x1_F1.iloc[0]
x_all[1]=x11_F1.iloc[0]
x_all[2]=x12_F1.iloc[0]
x_all[3]=x13_F1.iloc[0]
x_all[4]=x14_F1.iloc[0]
x_disc_F1=x_all.mean()
y_all=np.zeros(5)
y_all[0]=y1_F1.iloc[0]
y_all[1]=y11_F1.iloc[0]
y_all[2]=y12_F1.iloc[0]
y_all[3]=y13_F1.iloc[0]
y_all[4]=y14_F1.iloc[0]
y_disc_F1=y_all.mean()

###the last point of curves in F1
y1_tip_F1=y1_F1.iloc[-1]
x1_tip_F1=x1_F1.iloc[-1]
y2_tip_F1=y11_F1.iloc[-1]
x2_tip_F1=x11_F1.iloc[-1]
y3_tip_F1=y12_F1.iloc[-1]
x3_tip_F1=x12_F1.iloc[-1]
y4_tip_F1=y13_F1.iloc[-1]
x4_tip_F1=x13_F1.iloc[-1]
y5_tip_F1=y14_F1.iloc[-1]
x5_tip_F1=x14_F1.iloc[-1]

###the first point of curves in F1
y1_base_F1=y1_F1.iloc[0]
x1_base_F1=x1_F1.iloc[0]
y2_base_F1=y11_F1.iloc[0]
x2_base_F1=x11_F1.iloc[0]
y3_base_F1=y12_F1.iloc[0]
x3_base_F1=x12_F1.iloc[0]
y4_base_F1=y13_F1.iloc[0]
x4_base_F1=x13_F1.iloc[0]
y5_base_F1=y14_F1.iloc[0]
x5_base_F1=x14_F1.iloc[0]


##disc new= disc points-disc points

x_disc_F1_n=x_disc_F1-x_disc_F1
y_disc_F1_n=y_disc_F1-y_disc_F1

(x_disc_F1_n, y_disc_F1_n)

## basal point new= basal point - disc point

y1_base_F1_n=y1_F1.iloc[0]-y_disc_F1
x1_base_F1_n=x1_F1.iloc[0]-x_disc_F1
y2_base_F1_n=y11_F1.iloc[0]-y_disc_F1
x2_base_F1_n=x11_F1.iloc[0]-x_disc_F1
y3_base_F1_n=y12_F1.iloc[0]-y_disc_F1
x3_base_F1_n=x12_F1.iloc[0]-x_disc_F1
y4_base_F1_n=y13_F1.iloc[0]-y_disc_F1
x4_base_F1_n=x13_F1.iloc[0]-x_disc_F1
y5_base_F1_n=y14_F1.iloc[0]-y_disc_F1
x5_base_F1_n=x14_F1.iloc[0]-x_disc_F1

x_base=[x1_base_F1_n, x2_base_F1_n, x3_base_F1_n, x4_base_F1_n, x5_base_F1_n]
y_base=[y1_base_F1_n, y2_base_F1_n, y3_base_F1_n, y4_base_F1_n, y5_base_F1_n]

##tip point new= tip point - disc point

y1_tip_F1_n=y1_F1.iloc[-1]-y_disc_F1
x1_tip_F1_n=x1_F1.iloc[-1]-x_disc_F1
y2_tip_F1_n=y11_F1.iloc[-1]-y_disc_F1
x2_tip_F1_n=x11_F1.iloc[-1]-x_disc_F1
y3_tip_F1_n=y12_F1.iloc[-1]-y_disc_F1
x3_tip_F1_n=x12_F1.iloc[-1]-x_disc_F1
y4_tip_F1_n=y13_F1.iloc[-1]-y_disc_F1
x4_tip_F1_n=x13_F1.iloc[-1]-x_disc_F1
y5_tip_F1_n=y14_F1.iloc[-1]-y_disc_F1
x5_tip_F1_n=x14_F1.iloc[-1]-x_disc_F1

x_tip=[x1_tip_F1_n, x2_tip_F1_n, x3_tip_F1_n, x4_tip_F1_n, x5_tip_F1_n]
y_tip=[y1_tip_F1_n, y2_tip_F1_n, y3_tip_F1_n, y4_tip_F1_n, y5_tip_F1_n]

## alpha is the angle of disc center to basal point

m_alpha1= y1_base_F1_n/x1_base_F1_n
m_alpha2= y2_base_F1_n/x2_base_F1_n
m_alpha3= y3_base_F1_n/x3_base_F1_n
m_alpha4= y4_base_F1_n/x4_base_F1_n
m_alpha5= y5_base_F1_n/x5_base_F1_n

alpha1= np.degrees(np.arctan(m_alpha1))
alpha2= np.degrees(np.arctan(m_alpha2))
alpha3= np.degrees(np.arctan(m_alpha3))
alpha4= np.degrees(np.arctan(m_alpha4))
alpha5= np.degrees(np.arctan(m_alpha5))

alpha_total=[alpha1, alpha2, alpha3, alpha4, alpha5]


## betha is the angle of disc center to tip point

m_betha1= y1_tip_F1_n/x1_tip_F1_n
m_betha2= y2_tip_F1_n/x2_tip_F1_n
m_betha3= y3_tip_F1_n/x3_tip_F1_n
m_betha4= y4_tip_F1_n/x4_tip_F1_n
m_betha5= y5_tip_F1_n/x5_tip_F1_n

betha1= np.degrees(np.arctan(m_betha1))
betha2= np.degrees(np.arctan(m_betha2))
betha3= np.degrees(np.arctan(m_betha3))
betha4= np.degrees(np.arctan(m_betha4))
betha5= np.degrees(np.arctan(m_betha5))

betha_total= [betha1, betha2, betha3, betha4, betha5]
#########################################################

alpha_total_n=np.zeros(5)

for j in range (5):
    if x_base[j]<0 and y_base[j]>=0:
        alpha_total_n[j]= 180+alpha_total[j]
    
    elif x_base[j]<0 and y_base[j]<0:
        alpha_total_n[j]= 180-alpha_total[j]
    
    elif  x_base[j]>0 and y_base[j]<0:
        alpha_total_n[j]= -alpha_total[j]
        
    elif x_base[j]>0 and y_base[j]>0:
        alpha_total_n[j]=alpha_total[j]

print (alpha_total_n)

#######################################################
betha_total_n=np.zeros(5)
   
for k in range (5):
    if x_tip[k]<0 and y_tip[k]>=0:
        betha_total_n[k]= 180+betha_total[k]
    
    elif x_tip[k]<0 and y_tip[k]<0:
        betha_total_n[k]= 180-betha_total[k]
    
    elif x_tip[k]>0 and y_tip[k]<0:
        betha_total_n[k]= -betha_total[k]
        
    elif x_tip[k]>0 and y_tip[k]>0:
        betha_total_n[k]=betha_total[k]
        
print(betha_total_n)

######################################################

#### the difference between two lines

phi1= betha_total_n[0]-alpha_total_n[0]
phi2= betha_total_n[1]-alpha_total_n[1]
phi3= betha_total_n[2]-alpha_total_n[2]
phi4= betha_total_n[3]-alpha_total_n[3]
phi5= betha_total_n[4]-alpha_total_n[4]

phi_total_F1= [phi1, phi2, phi3, phi4, phi5]
phi_total_F1
##################################################

for i in range(5):
    if phi_total_F1[i] > 0:
        if phi_total_F1[i] < 180:
            print ('phi positive') 
        else:
            print('phi negative') 
    else:
        if abs(phi_total_F1[i]) <180:
            print('phi negative') 
        else:
            print('phi positive')    
            
##### ASSUMPTION: ABS PHI <180


# ### Frame 2

# In[ ]:



######choosing the rows for the first arm########
x1_F2=df_2['X-Coordinate '] [df_2['Curve Name']=='CURVE 1']
y1_F2=df_2['Y-Coordinate '] [df_2['Curve Name']=='CURVE 1']

######axes transformation########
x121=np.array((x1_F2-float(x1_F2[:1])))
y121=np.array((y1_F2-float(y1_F2[:1])))

######axes rotation########
m=(y121[-1])/x121[-1]
angle1 = np.degrees(np.arctan(m))

rules1 = [x121[-1]>0, y121[-1]>0, m>0]
rules2 = [x121[-1]>0, y121[-1]<0, m<0]
rules3 = [x121[-1]<0, y121[-1]<0, m>0]
rules4 = [x121[-1]<0, y121[-1]>0, m<0]

if  all(rules1):
    theta1= 360-angle1
elif  all(rules2): 
    theta1= abs(angle1)
elif  all(rules3): 
    theta1= 180-angle1
elif  all(rules4): 
    theta1= 180-angle1

x_F2=x121*(math.cos(math.radians(theta1)))-y121*(math.sin(math.radians(theta1)))
y_F2=x121*(math.sin(math.radians(theta1)))+y121*(math.cos(math.radians(theta1)))

######equidistant landmarks########
f=interp1d(x_F2,y_F2)
distance = np.cumsum(np.sqrt( np.ediff1d(x_F2, to_begin=0)**2 + np.ediff1d(y_F2, to_begin=0)**2 ))
distance = distance/distance[-1]

fx, fy = interp1d( distance, x_F2 ), interp1d( distance, y_F2 )

alpha = np.linspace(0, 1, 15)
x_regular_F2, y_regular_F2 = fx(alpha), fy(alpha)


################################to calculate sinuosity##################################

#calculating the length of a curve
d=np.zeros(14)
for i in range(14):
    d[i]=np.sqrt((x_regular_F2[i+1]-x_regular_F2[i])**2+(y_regular_F2[i+1]-y_regular_F2[i])**2)
print(d)
length=sum(d)

#calculating the distance between first and the last point (shotest length)
distance=np.sqrt((x_regular_F2[14]-x_regular_F2[0])**2+(y_regular_F2[14]-y_regular_F2[0])**2)

#calculating sinuosity
sinuosity1_F2=length/distance

#------------------------------------
x11_F2=df_2['X-Coordinate '] [df_2['Curve Name']=='CURVE 2']
y11_F2=df_2['Y-Coordinate '] [df_2['Curve Name']=='CURVE 2']
x122=np.array((x11_F2-float(x11_F2[:1])))
y122=np.array((y11_F2-float(y11_F2[:1])))


m2=(y122[-1])/x122[-1]
angle2 = np.degrees(np.arctan(m2))

rules1 = [x122[-1]>0, y122[-1]>0, m2>0]
rules2 = [x122[-1]>0, y122[-1]<0, m2<0]
rules3 = [x122[-1]<0, y122[-1]<0, m2>0]
rules4 = [x122[-1]<0, y122[-1]>0, m2<0]

if  all(rules1):
    theta2= 360-angle2
elif  all(rules2): 
    theta2= abs(angle2)
elif  all(rules3): 
    theta2= 180-angle2
elif  all(rules4): 
    theta2= 180-angle2
    
x2_F2=x122*(math.cos(math.radians(theta2)))-y122*(math.sin(math.radians(theta2)))
y2_F2=x122*(math.sin(math.radians(theta2)))+y122*(math.cos(math.radians(theta2)))

f=interp1d(x2_F2,y2_F2)

distance = np.cumsum(np.sqrt( np.ediff1d(x2_F2, to_begin=0)**2 + np.ediff1d(y2_F2, to_begin=0)**2 ))
distance = distance/distance[-1]

fx2, fy2 = interp1d( distance, x2_F2 ), interp1d( distance, y2_F2 )

alpha = np.linspace(0, 1, 15)
x_regular2_F2, y_regular2_F2 = fx2(alpha), fy2(alpha)


#to calculate sinuosity
d2=np.zeros(14)
for i in range(14):
    d2[i]=np.sqrt((x_regular2_F2[i+1]-x_regular2_F2[i])**2+(y_regular2_F2[i+1]-y_regular2_F2[i])**2)
print(d2)

length2=sum(d2)
distance2=np.sqrt((x_regular2_F2[14]-x_regular2_F2[0])**2+(y_regular2_F2[14]-y_regular2_F2[0])**2)
sinuosity2_F2=length2/distance2
#---------------------------------------------------------
x12_F2=df_2['X-Coordinate '] [df_2['Curve Name']=='CURVE 3']
y12_F2=df_2['Y-Coordinate '] [df_2['Curve Name']=='CURVE 3']
x123=np.array((x12_F2-float(x12_F2[:1])))
y123=np.array((y12_F2-float(y12_F2[:1])))

m3=(y123[-1])/x123[-1]
angle3 = np.degrees(np.arctan(m3))

rules1 = [x123[-1]>0, y123[-1]>0, m3>0]
rules2 = [x123[-1]>0, y123[-1]<0, m3<0]
rules3 = [x123[-1]<0, y123[-1]<0, m3>0]
rules4 = [x123[-1]<0, y123[-1]>0, m3<0]


if  all(rules1):
    theta3= 360-angle3
elif  all(rules2): 
    theta3= abs(angle3)
elif  all(rules3): 
    theta3= 180-angle3
elif  all(rules4): 
    theta3= 180-angle3


x3_F2=x123*(math.cos(math.radians(theta3)))-y123*(math.sin(math.radians(theta3)))
y3_F2=x123*(math.sin(math.radians(theta3)))+y123*(math.cos(math.radians(theta3)))


f=interp1d(x3_F2,y3_F2)
distance = np.cumsum(np.sqrt( np.ediff1d(x3_F2, to_begin=0)**2 + np.ediff1d(y3_F2, to_begin=0)**2 ))
distance = distance/distance[-1]

fx3, fy3 = interp1d( distance, x3_F2 ), interp1d( distance, y3_F2 )

alpha = np.linspace(0, 1, 15)
x_regular3_F2, y_regular3_F2 = fx3(alpha), fy3(alpha)

#to calculate sinuosity
d3=np.zeros(14)
for i in range(14):
    d3[i]=np.sqrt((x_regular3_F2[i+1]-x_regular3_F2[i])**2+(y_regular3_F2[i+1]-y_regular3_F2[i])**2)
print(d3)

length3=sum(d3)
distance3=np.sqrt((x_regular3_F2[14]-x_regular3_F2[0])**2+(y_regular3_F2[14]-y_regular3_F2[0])**2)
sinuosity3_F2=length3/distance3
#--------------------------------------------------
x13_F2=df_2['X-Coordinate '] [df_2['Curve Name']=='CURVE 4']
y13_F2=df_2['Y-Coordinate '] [df_2['Curve Name']=='CURVE 4']
x124=np.array((x13_F2-float(x13_F2[:1])))
y124=np.array((y13_F2-float(y13_F2[:1])))

m4=(y124[-1])/x124[-1]
angle4 = np.degrees(np.arctan(m4))

rules1 = [x124[-1]>0, y124[-1]>0, m4>0]
rules2 = [x124[-1]>0, y124[-1]<0, m4<0]
rules3 = [x124[-1]<0, y124[-1]<0, m4>0]
rules4 = [x124[-1]<0, y124[-1]>0, m4<0]


if  all(rules1):
    theta4= 360-angle4
elif  all(rules2): 
    theta4= abs(angle4)
elif  all(rules3): 
    theta4= 180-angle4
elif  all(rules4): 
    theta4= 180-angle4


x4_F2=x124*(math.cos(math.radians(theta4)))-y124*(math.sin(math.radians(theta4)))
y4_F2=x124*(math.sin(math.radians(theta4)))+y124*(math.cos(math.radians(theta4)))


f=interp1d(x4_F2,y4_F2)
distance = np.cumsum(np.sqrt( np.ediff1d(x4_F2, to_begin=0)**2 + np.ediff1d(y4_F2, to_begin=0)**2 ))
distance = distance/distance[-1]

fx4, fy4 = interp1d( distance, x4_F2 ), interp1d( distance, y4_F2 )

alpha = np.linspace(0, 1, 15)
x_regular4_F2, y_regular4_F2 = fx4(alpha), fy4(alpha)


#to calculate sinuosity
d4=np.zeros(15)
for i in range(14):
    d4[i]=np.sqrt((x_regular4_F2[i+1]-x_regular4_F2[i])**2+(y_regular4_F2[i+1]-y_regular4_F2[i])**2)
print(d4)

length4=sum(d4)
distance4=np.sqrt((x_regular4_F2[14]-x_regular4_F2[0])**2+(y_regular4_F2[14]-y_regular4_F2[0])**2)
sinuosity4_F2=length4/distance4
#--------------------------------------------------
x14_F2=df_2['X-Coordinate '] [df_2['Curve Name']=='CURVE 5']
y14_F2=df_2['Y-Coordinate '] [df_2['Curve Name']=='CURVE 5']
x125=np.array((x14_F2-float(x14_F2[:1])))
y125=np.array((y14_F2-float(y14_F2[:1])))


m5=(y125[-1])/x125[-1]
angle5= np.degrees(np.arctan(m5))

rules1 = [x125[-1]>0, y125[-1]>0, m5>0]
rules2 = [x125[-1]>0, y125[-1]<0, m5<0]
rules3 = [x125[-1]<0, y125[-1]<0, m5>0]
rules4 = [x125[-1]<0, y125[-1]>0, m5<0]


if  all(rules1):
    theta5= 360-angle5
elif  all(rules2): 
    theta5= abs(angle5)
elif  all(rules3): 
    theta5= 180-angle5
elif  all(rules4): 
    theta5= 180-angle5

x5_F2=x125*(math.cos(math.radians(theta5)))-y125*(math.sin(math.radians(theta5)))
y5_F2=x125*(math.sin(math.radians(theta5)))+y125*(math.cos(math.radians(theta5)))

f=interp1d(x5_F2,y5_F2)
distance = np.cumsum(np.sqrt( np.ediff1d(x5_F2, to_begin=0)**2 + np.ediff1d(y5_F2, to_begin=0)**2 ))
distance = distance/distance[-1]

fx5, fy5 = interp1d( distance, x5_F2), interp1d( distance, y5_F2 )

alpha = np.linspace(0, 1, 15)
x_regular5_F2, y_regular5_F2 = fx5(alpha), fy5(alpha)

#to calculate sinuosity
d5=np.zeros(14)
for i in range(14):
    d5[i]=np.sqrt((x_regular5_F2[i+1]-x_regular5_F2[i])**2+(y_regular5_F2[i+1]-y_regular5_F2[i])**2)
print(d5)

length5=sum(d5)
distance5=np.sqrt((x_regular5_F2[14]-x_regular5_F2[0])**2+(y_regular5_F2[14]-y_regular5_F2[0])**2)
sinuosity5_F2=length5/distance5

#center point##################################
x_all=np.zeros(5)
x_all[0]=x1_F2.iloc[0]
x_all[1]=x11_F2.iloc[0]
x_all[2]=x12_F2.iloc[0]
x_all[3]=x13_F2.iloc[0]
x_all[4]=x14_F2.iloc[0]
x_disc_F2=x_all.mean()
y_all=np.zeros(5)
y_all[0]=y1_F2.iloc[0]
y_all[1]=y11_F2.iloc[0]
y_all[2]=y12_F2.iloc[0]
y_all[3]=y13_F2.iloc[0]
y_all[4]=y14_F2.iloc[0]
y_disc_F2=y_all.mean()


###the last point of curves in  F2
y1_tip_F2=y1_F2.iloc[-1]
x1_tip_F2=x1_F2.iloc[-1]
y2_tip_F2=y11_F2.iloc[-1]
x2_tip_F2=x11_F2.iloc[-1]
y3_tip_F2=y12_F2.iloc[-1]
x3_tip_F2=x12_F2.iloc[-1]
y4_tip_F2=y13_F2.iloc[-1]
x4_tip_F2=x13_F2.iloc[-1]
y5_tip_F2=y14_F2.iloc[-1]
x5_tip_F2=x14_F2.iloc[-1]

###the first point of curves in F2
y1_base_F2=y1_F2.iloc[0]
x1_base_F2=x1_F2.iloc[0]
y2_base_F2=y11_F2.iloc[0]
x2_base_F2=x11_F2.iloc[0]
y3_base_F2=y12_F2.iloc[0]
x3_base_F2=x12_F2.iloc[0]
y4_base_F2=y13_F2.iloc[0]
x4_base_F2=x13_F2.iloc[0]
y5_base_F2=y14_F2.iloc[0]
x5_base_F2=x14_F2.iloc[0]


###################################################################### Arm angle F2

##disc new= disc points-disc points

x_disc_F2_n=x_disc_F2-x_disc_F2
y_disc_F2_n=y_disc_F2-y_disc_F2

(x_disc_F2_n, y_disc_F2_n)

## basal point new= basal point - disc point

y1_base_F2_n=y1_F2.iloc[0]-y_disc_F2
x1_base_F2_n=x1_F2.iloc[0]-x_disc_F2
y2_base_F2_n=y11_F2.iloc[0]-y_disc_F2
x2_base_F2_n=x11_F2.iloc[0]-x_disc_F2
y3_base_F2_n=y12_F2.iloc[0]-y_disc_F2
x3_base_F2_n=x12_F2.iloc[0]-x_disc_F2
y4_base_F2_n=y13_F2.iloc[0]-y_disc_F2
x4_base_F2_n=x13_F2.iloc[0]-x_disc_F2
y5_base_F2_n=y14_F2.iloc[0]-y_disc_F2
x5_base_F2_n=x14_F2.iloc[0]-x_disc_F2

x_base=[x1_base_F2_n, x2_base_F2_n, x3_base_F2_n, x4_base_F2_n, x5_base_F2_n]
y_base=[y1_base_F2_n, y2_base_F2_n, y3_base_F2_n, y4_base_F2_n, y5_base_F2_n]

##tip point new= tip point - disc point

y1_tip_F2_n=y1_F2.iloc[-1]-y_disc_F2
x1_tip_F2_n=x1_F2.iloc[-1]-x_disc_F2
y2_tip_F2_n=y11_F2.iloc[-1]-y_disc_F2
x2_tip_F2_n=x11_F2.iloc[-1]-x_disc_F2
y3_tip_F2_n=y12_F2.iloc[-1]-y_disc_F2
x3_tip_F2_n=x12_F2.iloc[-1]-x_disc_F2
y4_tip_F2_n=y13_F2.iloc[-1]-y_disc_F2
x4_tip_F2_n=x13_F2.iloc[-1]-x_disc_F2
y5_tip_F2_n=y14_F2.iloc[-1]-y_disc_F2
x5_tip_F2_n=x14_F2.iloc[-1]-x_disc_F2

x_tip=[x1_tip_F2_n, x2_tip_F2_n, x3_tip_F2_n, x4_tip_F2_n, x5_tip_F2_n]
y_tip=[y1_tip_F2_n, y2_tip_F2_n, y3_tip_F2_n, y4_tip_F2_n, y5_tip_F2_n]

## alpha is the angle of disc center to basal point

m_alpha1= y1_base_F2_n/x1_base_F2_n
m_alpha2= y2_base_F2_n/x2_base_F2_n
m_alpha3= y3_base_F2_n/x3_base_F2_n
m_alpha4= y4_base_F2_n/x4_base_F2_n
m_alpha5= y5_base_F2_n/x5_base_F2_n

alpha1= np.degrees(np.arctan(m_alpha1))
alpha2= np.degrees(np.arctan(m_alpha2))
alpha3= np.degrees(np.arctan(m_alpha3))
alpha4= np.degrees(np.arctan(m_alpha4))
alpha5= np.degrees(np.arctan(m_alpha5))

alpha_total=[alpha1, alpha2, alpha3, alpha4, alpha5]


## betha is the angle of disc center to tip point

m_betha1= y1_tip_F2_n/x1_tip_F2_n
m_betha2= y2_tip_F2_n/x2_tip_F2_n
m_betha3= y3_tip_F2_n/x3_tip_F2_n
m_betha4= y4_tip_F2_n/x4_tip_F2_n
m_betha5= y5_tip_F2_n/x5_tip_F2_n

betha1= np.degrees(np.arctan(m_betha1))
betha2= np.degrees(np.arctan(m_betha2))
betha3= np.degrees(np.arctan(m_betha3))
betha4= np.degrees(np.arctan(m_betha4))
betha5= np.degrees(np.arctan(m_betha5))

betha_total= [betha1, betha2, betha3, betha4, betha5]
#########################################################

alpha_total_n=np.zeros(5)

for j in range (5):
    if x_base[j]<0 and y_base[j]>0:
        alpha_total_n[j]= 180+alpha_total[j]
    
    elif x_base[j]<0 and y_base[j]<0:
        alpha_total_n[j]= 180-alpha_total[j]
    
    elif  x_base[j]>0 and y_base[j]<0:
        alpha_total_n[j]= -alpha_total[j]
        
    elif x_base[j]>0 and y_base[j]>0:
        alpha_total_n[j]=alpha_total[j]

print (alpha_total_n)
#######################################################
betha_total_n=np.zeros(5)
   
for k in range (5):
    if x_tip[k]<0 and y_tip[k]>0:
        betha_total_n[k]= 180+betha_total[k]
    
    elif x_tip[k]<0 and y_tip[k]<0:
        betha_total_n[k]= 180-betha_total[k]
    
    elif x_tip[k]>0 and y_tip[k]<0:
        betha_total_n[k]= -betha_total[k]
        
    elif x_tip[k]>0 and y_tip[k]>0:
        betha_total_n[k]=betha_total[k]
        
print(betha_total_n)

######################################################

#### the difference between two lines

phi1= betha_total_n[0]-alpha_total_n[0]
phi2= betha_total_n[1]-alpha_total_n[1]
phi3= betha_total_n[2]-alpha_total_n[2]
phi4= betha_total_n[3]-alpha_total_n[3]
phi5= betha_total_n[4]-alpha_total_n[4]

phi_total_F2= [phi1, phi2, phi3, phi4, phi5]
phi_total_F2

##################################################

for i in range(5):
    if phi_total_F2[i] > 0:
        if phi_total_F2[i] < 180:
            print ('phi positive') 
        else:
            print('phi negative') 
    else:
        if abs(phi_total_F2[i]) <180:
            print('phi negative') 
        else:
            print('phi positive')    
            
##### ASSUMPTION: ABS PHI <180


# ### Frame 3

# In[ ]:



######choosing the rows for the first arm########
x1_F2=df_2['X-Coordinate '] [df_2['Curve Name']=='CURVE 1']
y1_F2=df_2['Y-Coordinate '] [df_2['Curve Name']=='CURVE 1']

######axes transformation########
x121=np.array((x1_F2-float(x1_F2[:1])))
y121=np.array((y1_F2-float(y1_F2[:1])))

######axes rotation########
m=(y121[-1])/x121[-1]
angle1 = np.degrees(np.arctan(m))

rules1 = [x121[-1]>0, y121[-1]>0, m>0]
rules2 = [x121[-1]>0, y121[-1]<0, m<0]
rules3 = [x121[-1]<0, y121[-1]<0, m>0]
rules4 = [x121[-1]<0, y121[-1]>0, m<0]

if  all(rules1):
    theta1= 360-angle1
elif  all(rules2): 
    theta1= abs(angle1)
elif  all(rules3): 
    theta1= 180-angle1
elif  all(rules4): 
    theta1= 180-angle1

x_F2=x121*(math.cos(math.radians(theta1)))-y121*(math.sin(math.radians(theta1)))
y_F2=x121*(math.sin(math.radians(theta1)))+y121*(math.cos(math.radians(theta1)))

######eauidistant landmarks########
f=interp1d(x_F2,y_F2)
distance = np.cumsum(np.sqrt( np.ediff1d(x_F2, to_begin=0)**2 + np.ediff1d(y_F2, to_begin=0)**2 ))
distance = distance/distance[-1]

fx, fy = interp1d( distance, x_F2 ), interp1d( distance, y_F2 )

alpha = np.linspace(0, 1, 15)
x_regular_F2, y_regular_F2 = fx(alpha), fy(alpha)


################################to calculate sinuosity##################################

#calculating the length of a curve
d=np.zeros(14)
for i in range(14):
    d[i]=np.sqrt((x_regular_F2[i+1]-x_regular_F2[i])**2+(y_regular_F2[i+1]-y_regular_F2[i])**2)
print(d)
length=sum(d)

#calculating the distance between first and the last point (shotest length)
distance=np.sqrt((x_regular_F2[14]-x_regular_F2[0])**2+(y_regular_F2[14]-y_regular_F2[0])**2)

#calculating siniosity
sinuosity1_F2=length/distance

#------------------------------------
x11_F2=df_2['X-Coordinate '] [df_2['Curve Name']=='CURVE 2']
y11_F2=df_2['Y-Coordinate '] [df_2['Curve Name']=='CURVE 2']
x122=np.array((x11_F2-float(x11_F2[:1])))
y122=np.array((y11_F2-float(y11_F2[:1])))


m2=(y122[-1])/x122[-1]
angle2 = np.degrees(np.arctan(m2))

rules1 = [x122[-1]>0, y122[-1]>0, m2>0]
rules2 = [x122[-1]>0, y122[-1]<0, m2<0]
rules3 = [x122[-1]<0, y122[-1]<0, m2>0]
rules4 = [x122[-1]<0, y122[-1]>0, m2<0]

if  all(rules1):
    theta2= 360-angle2
elif  all(rules2): 
    theta2= abs(angle2)
elif  all(rules3): 
    theta2= 180-angle2
elif  all(rules4): 
    theta2= 180-angle2
    
x2_F2=x122*(math.cos(math.radians(theta2)))-y122*(math.sin(math.radians(theta2)))
y2_F2=x122*(math.sin(math.radians(theta2)))+y122*(math.cos(math.radians(theta2)))

f=interp1d(x2_F2,y2_F2)

distance = np.cumsum(np.sqrt( np.ediff1d(x2_F2, to_begin=0)**2 + np.ediff1d(y2_F2, to_begin=0)**2 ))
distance = distance/distance[-1]

fx2, fy2 = interp1d( distance, x2_F2 ), interp1d( distance, y2_F2 )

alpha = np.linspace(0, 1, 15)
x_regular2_F2, y_regular2_F2 = fx2(alpha), fy2(alpha)


#to calculate sinuosity
d2=np.zeros(14)
for i in range(14):
    d2[i]=np.sqrt((x_regular2_F2[i+1]-x_regular2_F2[i])**2+(y_regular2_F2[i+1]-y_regular2_F2[i])**2)
print(d2)

length2=sum(d2)
distance2=np.sqrt((x_regular2_F2[14]-x_regular2_F2[0])**2+(y_regular2_F2[14]-y_regular2_F2[0])**2)
sinuosity2_F2=length2/distance2
#---------------------------------------------------------
x12_F2=df_2['X-Coordinate '] [df_2['Curve Name']=='CURVE 3']
y12_F2=df_2['Y-Coordinate '] [df_2['Curve Name']=='CURVE 3']
x123=np.array((x12_F2-float(x12_F2[:1])))
y123=np.array((y12_F2-float(y12_F2[:1])))

m3=(y123[-1])/x123[-1]
angle3 = np.degrees(np.arctan(m3))

rules1 = [x123[-1]>0, y123[-1]>0, m3>0]
rules2 = [x123[-1]>0, y123[-1]<0, m3<0]
rules3 = [x123[-1]<0, y123[-1]<0, m3>0]
rules4 = [x123[-1]<0, y123[-1]>0, m3<0]


if  all(rules1):
    theta3= 360-angle3
elif  all(rules2): 
    theta3= abs(angle3)
elif  all(rules3): 
    theta3= 180-angle3
elif  all(rules4): 
    theta3= 180-angle3


x3_F2=x123*(math.cos(math.radians(theta3)))-y123*(math.sin(math.radians(theta3)))
y3_F2=x123*(math.sin(math.radians(theta3)))+y123*(math.cos(math.radians(theta3)))


f=interp1d(x3_F2,y3_F2)
distance = np.cumsum(np.sqrt( np.ediff1d(x3_F2, to_begin=0)**2 + np.ediff1d(y3_F2, to_begin=0)**2 ))
distance = distance/distance[-1]

fx3, fy3 = interp1d( distance, x3_F2 ), interp1d( distance, y3_F2 )

alpha = np.linspace(0, 1, 15)
x_regular3_F2, y_regular3_F2 = fx3(alpha), fy3(alpha)

#to calculate sinuosity
d3=np.zeros(14)
for i in range(14):
    d3[i]=np.sqrt((x_regular3_F2[i+1]-x_regular3_F2[i])**2+(y_regular3_F2[i+1]-y_regular3_F2[i])**2)
print(d3)

length3=sum(d3)
distance3=np.sqrt((x_regular3_F2[14]-x_regular3_F2[0])**2+(y_regular3_F2[14]-y_regular3_F2[0])**2)
sinuosity3_F2=length3/distance3
#--------------------------------------------------
x13_F2=df_2['X-Coordinate '] [df_2['Curve Name']=='CURVE 4']
y13_F2=df_2['Y-Coordinate '] [df_2['Curve Name']=='CURVE 4']
x124=np.array((x13_F2-float(x13_F2[:1])))
y124=np.array((y13_F2-float(y13_F2[:1])))

m4=(y124[-1])/x124[-1]
angle4 = np.degrees(np.arctan(m4))

rules1 = [x124[-1]>0, y124[-1]>0, m4>0]
rules2 = [x124[-1]>0, y124[-1]<0, m4<0]
rules3 = [x124[-1]<0, y124[-1]<0, m4>0]
rules4 = [x124[-1]<0, y124[-1]>0, m4<0]


if  all(rules1):
    theta4= 360-angle4
elif  all(rules2): 
    theta4= abs(angle4)
elif  all(rules3): 
    theta4= 180-angle4
elif  all(rules4): 
    theta4= 180-angle4


x4_F2=x124*(math.cos(math.radians(theta4)))-y124*(math.sin(math.radians(theta4)))
y4_F2=x124*(math.sin(math.radians(theta4)))+y124*(math.cos(math.radians(theta4)))


f=interp1d(x4_F2,y4_F2)
distance = np.cumsum(np.sqrt( np.ediff1d(x4_F2, to_begin=0)**2 + np.ediff1d(y4_F2, to_begin=0)**2 ))
distance = distance/distance[-1]

fx4, fy4 = interp1d( distance, x4_F2 ), interp1d( distance, y4_F2 )

alpha = np.linspace(0, 1, 15)
x_regular4_F2, y_regular4_F2 = fx4(alpha), fy4(alpha)


#to calculate sinuosity
d4=np.zeros(15)
for i in range(14):
    d4[i]=np.sqrt((x_regular4_F2[i+1]-x_regular4_F2[i])**2+(y_regular4_F2[i+1]-y_regular4_F2[i])**2)
print(d4)

length4=sum(d4)
distance4=np.sqrt((x_regular4_F2[14]-x_regular4_F2[0])**2+(y_regular4_F2[14]-y_regular4_F2[0])**2)
sinuosity4_F2=length4/distance4
#--------------------------------------------------
x14_F2=df_2['X-Coordinate '] [df_2['Curve Name']=='CURVE 5']
y14_F2=df_2['Y-Coordinate '] [df_2['Curve Name']=='CURVE 5']
x125=np.array((x14_F2-float(x14_F2[:1])))
y125=np.array((y14_F2-float(y14_F2[:1])))


m5=(y125[-1])/x125[-1]
angle5= np.degrees(np.arctan(m5))

rules1 = [x125[-1]>0, y125[-1]>0, m5>0]
rules2 = [x125[-1]>0, y125[-1]<0, m5<0]
rules3 = [x125[-1]<0, y125[-1]<0, m5>0]
rules4 = [x125[-1]<0, y125[-1]>0, m5<0]


if  all(rules1):
    theta5= 360-angle5
elif  all(rules2): 
    theta5= abs(angle5)
elif  all(rules3): 
    theta5= 180-angle5
elif  all(rules4): 
    theta5= 180-angle5

x5_F2=x125*(math.cos(math.radians(theta5)))-y125*(math.sin(math.radians(theta5)))
y5_F2=x125*(math.sin(math.radians(theta5)))+y125*(math.cos(math.radians(theta5)))

f=interp1d(x5_F2,y5_F2)
distance = np.cumsum(np.sqrt( np.ediff1d(x5_F2, to_begin=0)**2 + np.ediff1d(y5_F2, to_begin=0)**2 ))
distance = distance/distance[-1]

fx5, fy5 = interp1d( distance, x5_F2), interp1d( distance, y5_F2 )

alpha = np.linspace(0, 1, 15)
x_regular5_F2, y_regular5_F2 = fx5(alpha), fy5(alpha)

#to calculate sinuosity
d5=np.zeros(14)
for i in range(14):
    d5[i]=np.sqrt((x_regular5_F2[i+1]-x_regular5_F2[i])**2+(y_regular5_F2[i+1]-y_regular5_F2[i])**2)
print(d5)

length5=sum(d5)
distance5=np.sqrt((x_regular5_F2[14]-x_regular5_F2[0])**2+(y_regular5_F2[14]-y_regular5_F2[0])**2)
sinuosity5_F2=length5/distance5

#center point##################################
x_all=np.zeros(5)
x_all[0]=x1_F2.iloc[0]
x_all[1]=x11_F2.iloc[0]
x_all[2]=x12_F2.iloc[0]
x_all[3]=x13_F2.iloc[0]
x_all[4]=x14_F2.iloc[0]
x_disc_F2=x_all.mean()
y_all=np.zeros(5)
y_all[0]=y1_F2.iloc[0]
y_all[1]=y11_F2.iloc[0]
y_all[2]=y12_F2.iloc[0]
y_all[3]=y13_F2.iloc[0]
y_all[4]=y14_F2.iloc[0]
y_disc_F2=y_all.mean()


###the last point of curves in  F2
y1_tip_F2=y1_F2.iloc[-1]
x1_tip_F2=x1_F2.iloc[-1]
y2_tip_F2=y11_F2.iloc[-1]
x2_tip_F2=x11_F2.iloc[-1]
y3_tip_F2=y12_F2.iloc[-1]
x3_tip_F2=x12_F2.iloc[-1]
y4_tip_F2=y13_F2.iloc[-1]
x4_tip_F2=x13_F2.iloc[-1]
y5_tip_F2=y14_F2.iloc[-1]
x5_tip_F2=x14_F2.iloc[-1]

###the first point of curves in F2
y1_base_F2=y1_F2.iloc[0]
x1_base_F2=x1_F2.iloc[0]
y2_base_F2=y11_F2.iloc[0]
x2_base_F2=x11_F2.iloc[0]
y3_base_F2=y12_F2.iloc[0]
x3_base_F2=x12_F2.iloc[0]
y4_base_F2=y13_F2.iloc[0]
x4_base_F2=x13_F2.iloc[0]
y5_base_F2=y14_F2.iloc[0]
x5_base_F2=x14_F2.iloc[0]


###################################################################### Arm angle F2

##disc new= disc points-disc points

x_disc_F2_n=x_disc_F2-x_disc_F2
y_disc_F2_n=y_disc_F2-y_disc_F2

(x_disc_F2_n, y_disc_F2_n)

## basal point new= basal point - disc point

y1_base_F2_n=y1_F2.iloc[0]-y_disc_F2
x1_base_F2_n=x1_F2.iloc[0]-x_disc_F2
y2_base_F2_n=y11_F2.iloc[0]-y_disc_F2
x2_base_F2_n=x11_F2.iloc[0]-x_disc_F2
y3_base_F2_n=y12_F2.iloc[0]-y_disc_F2
x3_base_F2_n=x12_F2.iloc[0]-x_disc_F2
y4_base_F2_n=y13_F2.iloc[0]-y_disc_F2
x4_base_F2_n=x13_F2.iloc[0]-x_disc_F2
y5_base_F2_n=y14_F2.iloc[0]-y_disc_F2
x5_base_F2_n=x14_F2.iloc[0]-x_disc_F2

x_base=[x1_base_F2_n, x2_base_F2_n, x3_base_F2_n, x4_base_F2_n, x5_base_F2_n]
y_base=[y1_base_F2_n, y2_base_F2_n, y3_base_F2_n, y4_base_F2_n, y5_base_F2_n]

##tip point new= tip point - disc point

y1_tip_F2_n=y1_F2.iloc[-1]-y_disc_F2
x1_tip_F2_n=x1_F2.iloc[-1]-x_disc_F2
y2_tip_F2_n=y11_F2.iloc[-1]-y_disc_F2
x2_tip_F2_n=x11_F2.iloc[-1]-x_disc_F2
y3_tip_F2_n=y12_F2.iloc[-1]-y_disc_F2
x3_tip_F2_n=x12_F2.iloc[-1]-x_disc_F2
y4_tip_F2_n=y13_F2.iloc[-1]-y_disc_F2
x4_tip_F2_n=x13_F2.iloc[-1]-x_disc_F2
y5_tip_F2_n=y14_F2.iloc[-1]-y_disc_F2
x5_tip_F2_n=x14_F2.iloc[-1]-x_disc_F2

x_tip=[x1_tip_F2_n, x2_tip_F2_n, x3_tip_F2_n, x4_tip_F2_n, x5_tip_F2_n]
y_tip=[y1_tip_F2_n, y2_tip_F2_n, y3_tip_F2_n, y4_tip_F2_n, y5_tip_F2_n]

## alpha is the angle of disc center to basal point

m_alpha1= y1_base_F2_n/x1_base_F2_n
m_alpha2= y2_base_F2_n/x2_base_F2_n
m_alpha3= y3_base_F2_n/x3_base_F2_n
m_alpha4= y4_base_F2_n/x4_base_F2_n
m_alpha5= y5_base_F2_n/x5_base_F2_n

alpha1= np.degrees(np.arctan(m_alpha1))
alpha2= np.degrees(np.arctan(m_alpha2))
alpha3= np.degrees(np.arctan(m_alpha3))
alpha4= np.degrees(np.arctan(m_alpha4))
alpha5= np.degrees(np.arctan(m_alpha5))

alpha_total=[alpha1, alpha2, alpha3, alpha4, alpha5]


## betha is the angle of disc center to tip point

m_betha1= y1_tip_F2_n/x1_tip_F2_n
m_betha2= y2_tip_F2_n/x2_tip_F2_n
m_betha3= y3_tip_F2_n/x3_tip_F2_n
m_betha4= y4_tip_F2_n/x4_tip_F2_n
m_betha5= y5_tip_F2_n/x5_tip_F2_n

betha1= np.degrees(np.arctan(m_betha1))
betha2= np.degrees(np.arctan(m_betha2))
betha3= np.degrees(np.arctan(m_betha3))
betha4= np.degrees(np.arctan(m_betha4))
betha5= np.degrees(np.arctan(m_betha5))

betha_total= [betha1, betha2, betha3, betha4, betha5]
#########################################################

alpha_total_n=np.zeros(5)

for j in range (5):
    if x_base[j]<0 and y_base[j]>0:
        alpha_total_n[j]= 180+alpha_total[j]
    
    elif x_base[j]<0 and y_base[j]<0:
        alpha_total_n[j]= 180-alpha_total[j]
    
    elif  x_base[j]>0 and y_base[j]<0:
        alpha_total_n[j]= -alpha_total[j]
        
    elif x_base[j]>0 and y_base[j]>0:
        alpha_total_n[j]=alpha_total[j]

print (alpha_total_n)
#######################################################
betha_total_n=np.zeros(5)
   
for k in range (5):
    if x_tip[k]<0 and y_tip[k]>0:
        betha_total_n[k]= 180+betha_total[k]
    
    elif x_tip[k]<0 and y_tip[k]<0:
        betha_total_n[k]= 180-betha_total[k]
    
    elif x_tip[k]>0 and y_tip[k]<0:
        betha_total_n[k]= -betha_total[k]
        
    elif x_tip[k]>0 and y_tip[k]>0:
        betha_total_n[k]=betha_total[k]
        
print(betha_total_n)

######################################################

#### the difference between two lines

phi1= betha_total_n[0]-alpha_total_n[0]
phi2= betha_total_n[1]-alpha_total_n[1]
phi3= betha_total_n[2]-alpha_total_n[2]
phi4= betha_total_n[3]-alpha_total_n[3]
phi5= betha_total_n[4]-alpha_total_n[4]

phi_total_F2= [phi1, phi2, phi3, phi4, phi5]
phi_total_F2

##################################################

for i in range(5):
    if phi_total_F2[i] > 0:
        if phi_total_F2[i] < 180:
            print ('phi positive') 
        else:
            print('phi negative') 
    else:
        if abs(phi_total_F2[i]) <180:
            print('phi negative') 
        else:
            print('phi positive')    
            
##### ASSUMPTION: ABS PHI <180


# ### Frame 4

# In[ ]:



######choosing the rows for the first arm########
x1_F4=df_4['X-Coordinate '] [df_4['Curve Name']=='CURVE 1']
y1_F4=df_4['Y-Coordinate '] [df_4['Curve Name']=='CURVE 1']

######axes transformation########
x121=np.array((x1_F4-float(x1_F4[:1])))
y121=np.array((y1_F4-float(y1_F4[:1])))

######axes rotation########
m=(y121[-1])/x121[-1]
angle1 = np.degrees(np.arctan(m))

rules1 = [x121[-1]>0, y121[-1]>0, m>0]
rules2 = [x121[-1]>0, y121[-1]<0, m<0]
rules3 = [x121[-1]<0, y121[-1]<0, m>0]
rules4 = [x121[-1]<0, y121[-1]>0, m<0]


if  all(rules1):
    theta1= 360-angle1
elif  all(rules2): 
    theta1= abs(angle1)
elif  all(rules3): 
    theta1= 180-angle1
elif  all(rules4): 
    theta1= 180-angle1

x_F4=x121*(math.cos(math.radians(theta1)))-y121*(math.sin(math.radians(theta1)))
y_F4=x121*(math.sin(math.radians(theta1)))+y121*(math.cos(math.radians(theta1)))

######equidistant landmarks########
f=interp1d(x_F4,y_F4)
distance = np.cumsum(np.sqrt( np.ediff1d(x_F4, to_begin=0)**2 + np.ediff1d(y_F4, to_begin=0)**2 ))
distance = distance/distance[-1]

fx, fy = interp1d( distance, x_F4 ), interp1d( distance, y_F4 )

alpha = np.linspace(0, 1, 15)
x_regular_F4, y_regular_F4 = fx(alpha), fy(alpha)


################################to calculate sinuosity##################################

#calculating the length of a curve
d=np.zeros(14)
for i in range(14):
    d[i]=np.sqrt((x_regular_F4[i+1]-x_regular_F4[i])**2+(y_regular_F4[i+1]-y_regular_F4[i])**2)
print(d)
length=sum(d)

#calculating the distance between first and the last point (shotest length)
distance=np.sqrt((x_regular_F4[14]-x_regular_F4[0])**2+(y_regular_F4[14]-y_regular_F4[0])**2)

#calculating sinuosity
sinuosity1_F4=length/distance

#------------------------------------
x11_F4=df_4['X-Coordinate '] [df_4['Curve Name']=='CURVE 2']
y11_F4=df_4['Y-Coordinate '] [df_4['Curve Name']=='CURVE 2']
x122=np.array((x11_F4-float(x11_F4[:1])))
y122=np.array((y11_F4-float(y11_F4[:1])))


m2=(y122[-1])/x122[-1]
angle2 = np.degrees(np.arctan(m2))

rules1 = [x122[-1]>0, y122[-1]>0, m2>0]
rules2 = [x122[-1]>0, y122[-1]<0, m2<0]
rules3 = [x122[-1]<0, y122[-1]<0, m2>0]
rules4 = [x122[-1]<0, y122[-1]>0, m2<0]

if  all(rules1):
    theta2= 360-angle2
elif  all(rules2): 
    theta2= abs(angle2)
elif  all(rules3): 
    theta2= 180-angle2
elif  all(rules4): 
    theta2= 180-angle2
    
x2_F4=x122*(math.cos(math.radians(theta2)))-y122*(math.sin(math.radians(theta2)))
y2_F4=x122*(math.sin(math.radians(theta2)))+y122*(math.cos(math.radians(theta2)))

f=interp1d(x2_F4,y2_F4)

distance = np.cumsum(np.sqrt( np.ediff1d(x2_F4, to_begin=0)**2 + np.ediff1d(y2_F4, to_begin=0)**2 ))
distance = distance/distance[-1]

fx2, fy2 = interp1d( distance, x2_F4 ), interp1d( distance, y2_F4 )

alpha = np.linspace(0, 1, 15)
x_regular2_F4, y_regular2_F4 = fx2(alpha), fy2(alpha)


#to calculate sinuosity
d2=np.zeros(14)
for i in range(14):
    d2[i]=np.sqrt((x_regular2_F4[i+1]-x_regular2_F4[i])**2+(y_regular2_F4[i+1]-y_regular2_F4[i])**2)
print(d2)

length2=sum(d2)
distance2=np.sqrt((x_regular2_F4[14]-x_regular2_F4[0])**2+(y_regular2_F4[14]-y_regular2_F4[0])**2)
sinuosity2_F4=length2/distance2
#---------------------------------------------------------
x12_F4=df_4['X-Coordinate '] [df_4['Curve Name']=='CURVE 3']
y12_F4=df_4['Y-Coordinate '] [df_4['Curve Name']=='CURVE 3']
x123=np.array((x12_F4-float(x12_F4[:1])))
y123=np.array((y12_F4-float(y12_F4[:1])))

m3=(y123[-1])/x123[-1]
angle3 = np.degrees(np.arctan(m3))

rules1 = [x123[-1]>0, y123[-1]>0, m3>0]
rules2 = [x123[-1]>0, y123[-1]<0, m3<0]
rules3 = [x123[-1]<0, y123[-1]<0, m3>0]
rules4 = [x123[-1]<0, y123[-1]>0, m3<0]


if  all(rules1):
    theta3= 360-angle3
elif  all(rules2): 
    theta3= abs(angle3)
elif  all(rules3): 
    theta3= 180-angle3
elif  all(rules4): 
    theta3= 180-angle3


x3_F4=x123*(math.cos(math.radians(theta3)))-y123*(math.sin(math.radians(theta3)))
y3_F4=x123*(math.sin(math.radians(theta3)))+y123*(math.cos(math.radians(theta3)))


f=interp1d(x3_F4,y3_F4)
distance = np.cumsum(np.sqrt( np.ediff1d(x3_F4, to_begin=0)**2 + np.ediff1d(y3_F4, to_begin=0)**2 ))
distance = distance/distance[-1]

fx3, fy3 = interp1d( distance, x3_F4 ), interp1d( distance, y3_F4 )

alpha = np.linspace(0, 1, 15)
x_regular3_F4, y_regular3_F4 = fx3(alpha), fy3(alpha)


#to calculate sinuosity
d3=np.zeros(14)
for i in range(14):
    d3[i]=np.sqrt((x_regular3_F4[i+1]-x_regular3_F4[i])**2+(y_regular3_F4[i+1]-y_regular3_F4[i])**2)
print(d3)

length3=sum(d3)
distance3=np.sqrt((x_regular3_F4[14]-x_regular3_F4[0])**2+(y_regular3_F4[14]-y_regular3_F4[0])**2)
sinuosity3_F4=length3/distance3
#--------------------------------------------------
x13_F4=df_4['X-Coordinate '] [df_4['Curve Name']=='CURVE 4']
y13_F4=df_4['Y-Coordinate '] [df_4['Curve Name']=='CURVE 4']
x124=np.array((x13_F4-float(x13_F4[:1])))
y124=np.array((y13_F4-float(y13_F4[:1])))

m4=(y124[-1])/x124[-1]
angle4 = np.degrees(np.arctan(m4))

rules1 = [x124[-1]>0, y124[-1]>0, m4>0]
rules2 = [x124[-1]>0, y124[-1]<0, m4<0]
rules3 = [x124[-1]<0, y124[-1]<0, m4>0]
rules4 = [x124[-1]<0, y124[-1]>0, m4<0]


if  all(rules1):
    theta4= 360-angle4
elif  all(rules2): 
    theta4= abs(angle4)
elif  all(rules3): 
    theta4= 180-angle4
elif  all(rules4): 
    theta4= 180-angle4


x4_F4=x124*(math.cos(math.radians(theta4)))-y124*(math.sin(math.radians(theta4)))
y4_F4=x124*(math.sin(math.radians(theta4)))+y124*(math.cos(math.radians(theta4)))


f=interp1d(x4_F4,y4_F4)
distance = np.cumsum(np.sqrt( np.ediff1d(x4_F4, to_begin=0)**2 + np.ediff1d(y4_F4, to_begin=0)**2 ))
distance = distance/distance[-1]

fx4, fy4 = interp1d( distance, x4_F4 ), interp1d( distance, y4_F4 )

alpha = np.linspace(0, 1, 15)
x_regular4_F4, y_regular4_F4 = fx4(alpha), fy4(alpha)


#to calculate sinuosity
d4=np.zeros(15)
for i in range(14):
    d4[i]=np.sqrt((x_regular4_F4[i+1]-x_regular4_F4[i])**2+(y_regular4_F4[i+1]-y_regular4_F4[i])**2)
print(d4)

length4=sum(d4)
distance4=np.sqrt((x_regular4_F4[14]-x_regular4_F4[0])**2+(y_regular4_F4[14]-y_regular4_F4[0])**2)
sinuosity4_F4=length4/distance4
#--------------------------------------------------
x14_F4=df_4['X-Coordinate '] [df_4['Curve Name']=='CURVE 5']
y14_F4=df_4['Y-Coordinate '] [df_4['Curve Name']=='CURVE 5']
x125=np.array((x14_F4-float(x14_F4[:1])))
y125=np.array((y14_F4-float(y14_F4[:1])))


m5=(y125[-1])/x125[-1]
angle5= np.degrees(np.arctan(m5))

rules1 = [x125[-1]>0, y125[-1]>0, m5>0]
rules2 = [x125[-1]>0, y125[-1]<0, m5<0]
rules3 = [x125[-1]<0, y125[-1]<0, m5>0]
rules4 = [x125[-1]<0, y125[-1]>0, m5<0]


if  all(rules1):
    theta5= 360-angle5
elif  all(rules2): 
    theta5= abs(angle5)
elif  all(rules3): 
    theta5= 180-angle5
elif  all(rules4): 
    theta5= 180-angle5


x5_F4=x125*(math.cos(math.radians(theta5)))-y125*(math.sin(math.radians(theta5)))
y5_F4=x125*(math.sin(math.radians(theta5)))+y125*(math.cos(math.radians(theta5)))

f=interp1d(x5_F4,y5_F4)
distance = np.cumsum(np.sqrt( np.ediff1d(x5_F4, to_begin=0)**2 + np.ediff1d(y5_F4, to_begin=0)**2 ))
distance = distance/distance[-1]

fx5, fy5 = interp1d( distance, x5_F4), interp1d( distance, y5_F4 )

alpha = np.linspace(0, 1, 15)
x_regular5_F4, y_regular5_F4 = fx5(alpha), fy5(alpha)

#to calculate sinuosity
d5=np.zeros(14)
for i in range(14):
    d5[i]=np.sqrt((x_regular5_F4[i+1]-x_regular5_F4[i])**2+(y_regular5_F4[i+1]-y_regular5_F4[i])**2)
print(d5)

length5=sum(d5)
distance5=np.sqrt((x_regular5_F4[14]-x_regular5_F4[0])**2+(y_regular5_F4[14]-y_regular5_F4[0])**2)
sinuosity5_F4=length5/distance5

#center point##################################
x_all=np.zeros(5)
x_all[0]=x1_F4.iloc[0]
x_all[1]=x11_F4.iloc[0]
x_all[2]=x12_F4.iloc[0]
x_all[3]=x13_F4.iloc[0]
x_all[4]=x14_F4.iloc[0]
x_disc_F4=x_all.mean()
y_all=np.zeros(5)
y_all[0]=y1_F4.iloc[0]
y_all[1]=y11_F4.iloc[0]
y_all[2]=y12_F4.iloc[0]
y_all[3]=y13_F4.iloc[0]
y_all[4]=y14_F4.iloc[0]
y_disc_F4=y_all.mean()


###the last point of curves in F4
y1_tip_F4=y1_F4.iloc[-1]
x1_tip_F4=x1_F4.iloc[-1]
y2_tip_F4=y11_F4.iloc[-1]
x2_tip_F4=x11_F4.iloc[-1]
y3_tip_F4=y12_F4.iloc[-1]
x3_tip_F4=x12_F4.iloc[-1]
y4_tip_F4=y13_F4.iloc[-1]
x4_tip_F4=x13_F4.iloc[-1]
y5_tip_F4=y14_F4.iloc[-1]
x5_tip_F4=x14_F4.iloc[-1]


###the first point of curves in F4
y1_base_F4=y1_F4.iloc[0]
x1_base_F4=x1_F4.iloc[0]
y2_base_F4=y11_F4.iloc[0]
x2_base_F4=x11_F4.iloc[0]
y3_base_F4=y12_F4.iloc[0]
x3_base_F4=x12_F4.iloc[0]
y4_base_F4=y13_F4.iloc[0]
x4_base_F4=x13_F4.iloc[0]
y5_base_F4=y14_F4.iloc[0]
x5_base_F4=x14_F4.iloc[0]


###################################################################### Arm angle F4

##disc new= disc points-disc points

x_disc_F4_n=x_disc_F4-x_disc_F4
y_disc_F4_n=y_disc_F4-y_disc_F4


## basal point new= basal point - disc point

y1_base_F4_n=y1_F4.iloc[0]-y_disc_F4
x1_base_F4_n=x1_F4.iloc[0]-x_disc_F4
y2_base_F4_n=y11_F4.iloc[0]-y_disc_F4
x2_base_F4_n=x11_F4.iloc[0]-x_disc_F4
y3_base_F4_n=y12_F4.iloc[0]-y_disc_F4
x3_base_F4_n=x12_F4.iloc[0]-x_disc_F4
y4_base_F4_n=y13_F4.iloc[0]-y_disc_F4
x4_base_F4_n=x13_F4.iloc[0]-x_disc_F4
y5_base_F4_n=y14_F4.iloc[0]-y_disc_F4
x5_base_F4_n=x14_F4.iloc[0]-x_disc_F4

x_base=[x1_base_F4_n, x2_base_F4_n, x3_base_F4_n, x4_base_F4_n, x5_base_F4_n]
y_base=[y1_base_F4_n, y2_base_F4_n, y3_base_F4_n, y4_base_F4_n, y5_base_F4_n]

##tip point new= tip point - disc point

y1_tip_F4_n=y1_F4.iloc[-1]-y_disc_F4
x1_tip_F4_n=x1_F4.iloc[-1]-x_disc_F4
y2_tip_F4_n=y11_F4.iloc[-1]-y_disc_F4
x2_tip_F4_n=x11_F4.iloc[-1]-x_disc_F4
y3_tip_F4_n=y12_F4.iloc[-1]-y_disc_F4
x3_tip_F4_n=x12_F4.iloc[-1]-x_disc_F4
y4_tip_F4_n=y13_F4.iloc[-1]-y_disc_F4
x4_tip_F4_n=x13_F4.iloc[-1]-x_disc_F4
y5_tip_F4_n=y14_F4.iloc[-1]-y_disc_F4
x5_tip_F4_n=x14_F4.iloc[-1]-x_disc_F4

x_tip=[x1_tip_F4_n, x2_tip_F4_n, x3_tip_F4_n, x4_tip_F4_n, x5_tip_F4_n]
y_tip=[y1_tip_F4_n, y2_tip_F4_n, y3_tip_F4_n, y4_tip_F4_n, y5_tip_F4_n]

## alpha is the angle of disc center to basal point

m_alpha1= y1_base_F4_n/x1_base_F4_n
m_alpha2= y2_base_F4_n/x2_base_F4_n
m_alpha3= y3_base_F4_n/x3_base_F4_n
m_alpha4= y4_base_F4_n/x4_base_F4_n
m_alpha5= y5_base_F4_n/x5_base_F4_n

alpha1= np.degrees(np.arctan(m_alpha1))
alpha2= np.degrees(np.arctan(m_alpha2))
alpha3= np.degrees(np.arctan(m_alpha3))
alpha4= np.degrees(np.arctan(m_alpha4))
alpha5= np.degrees(np.arctan(m_alpha5))

alpha_total=[alpha1, alpha2, alpha3, alpha4, alpha5]


## betha is the angle of disc center to tip point

m_betha1= y1_tip_F4_n/x1_tip_F4_n
m_betha2= y2_tip_F4_n/x2_tip_F4_n
m_betha3= y3_tip_F4_n/x3_tip_F4_n
m_betha4= y4_tip_F4_n/x4_tip_F4_n
m_betha5= y5_tip_F4_n/x5_tip_F4_n

betha1= np.degrees(np.arctan(m_betha1))
betha2= np.degrees(np.arctan(m_betha2))
betha3= np.degrees(np.arctan(m_betha3))
betha4= np.degrees(np.arctan(m_betha4))
betha5= np.degrees(np.arctan(m_betha5))

betha_total= [betha1, betha2, betha3, betha4, betha5]
#########################################################

alpha_total_n=np.zeros(5)

for j in range (5):
    if x_base[j]<0 and y_base[j]>0:
        alpha_total_n[j]= 180+alpha_total[j]
    
    elif x_base[j]<0 and y_base[j]<0:
        alpha_total_n[j]= 180-alpha_total[j]
    
    elif  x_base[j]>0 and y_base[j]<0:
        alpha_total_n[j]= -alpha_total[j]
        
    elif x_base[j]>0 and y_base[j]>0:
        alpha_total_n[j]=alpha_total[j]

print (alpha_total_n)
#######################################################
betha_total_n=np.zeros(5)
   
for k in range (5):
    if x_tip[k]<0 and y_tip[k]>0:
        betha_total_n[k]= 180+betha_total[k]
    
    elif x_tip[k]<0 and y_tip[k]<0:
        betha_total_n[k]= 180-betha_total[k]
    
    elif x_tip[k]>0 and y_tip[k]<0:
        betha_total_n[k]= -betha_total[k]
        
    elif x_tip[k]>0 and y_tip[k]>0:
        betha_total_n[k]=betha_total[k]
        
print(betha_total_n)

######################################################

#### the difference between two lines

phi1= betha_total_n[0]-alpha_total_n[0]
phi2= betha_total_n[1]-alpha_total_n[1]
phi3= betha_total_n[2]-alpha_total_n[2]
phi4= betha_total_n[3]-alpha_total_n[3]
phi5= betha_total_n[4]-alpha_total_n[4]

phi_total_F4= [phi1, phi2, phi3, phi4, phi5]
phi_total_F4

##################################################

for i in range(5):
    if phi_total_F4[i] > 0:
        if phi_total_F4[i] < 180:
            print ('phi positive') 
        else:
            print('phi negative') 
    else:
        if abs(phi_total_F4[i]) <180:
            print('phi negative') 
        else:
            print('phi positive')    
            
##### ASSUMPTION: ABS PHI <180


# ### Frame 5

# In[ ]:



######choosing the rows for the first arm########
x1_F5=df_5['X-Coordinate '] [df_5['Curve Name']=='CURVE 1']
y1_F5=(df_5['Y-Coordinate '] [df_5['Curve Name']=='CURVE 1'])

######axes transformation########
x121=np.array((x1_F5-float(x1_F5[:1])))
y121=np.array((y1_F5-float(y1_F5[:1])))

######axes rotation########
m=(y121[-1])/x121[-1]
angle1 = np.degrees(np.arctan(m))

rules1 = [x121[-1]>0, y121[-1]>0, m>0]
rules2 = [x121[-1]>0, y121[-1]<0, m<0]
rules3 = [x121[-1]<0, y121[-1]<0, m>0]
rules4 = [x121[-1]<0, y121[-1]>0, m<0]


if  all(rules1):
    theta1= 360-angle1
elif  all(rules2): 
    theta1= abs(angle1)
elif  all(rules3): 
    theta1= 180-angle1
elif  all(rules4): 
    theta1= 180-angle1
print(theta1)

x_F5=x121*(math.cos(math.radians(theta1)))-y121*(math.sin(math.radians(theta1)))
y_F5=x121*(math.sin(math.radians(theta1)))+y121*(math.cos(math.radians(theta1)))


######equidistant landmarks########
f=interp1d(x_F5,y_F5)
distance = np.cumsum(np.sqrt( np.ediff1d(x_F5, to_begin=0)**2 + np.ediff1d(y_F5, to_begin=0)**2 ))
distance = distance/distance[-1]

fx, fy = interp1d( distance, x_F5 ), interp1d( distance, y_F5 )

alpha = np.linspace(0, 1, 15)
x_regular_F5, y_regular_F5 = fx(alpha), fy(alpha)


################################to calculate sinuosity##################################

#calculating the length of a curve
d=np.zeros(14)
for i in range(14):
    d[i]=np.sqrt((x_regular_F5[i+1]-x_regular_F5[i])**2+(y_regular_F5[i+1]-y_regular_F5[i])**2)
print(d)
length=sum(d)

#calculating the distance between first and the last point (shotest length)
distance=np.sqrt((x_regular_F5[14]-x_regular_F5[0])**2+(y_regular_F5[14]-y_regular_F5[0])**2)

#calculating sinuosity
sinuosity1_F5=length/distance

#------------------------------------
x11_F5=df_5['X-Coordinate '] [df_5['Curve Name']=='CURVE 2']
y11_F5=(df_5['Y-Coordinate '] [df_5['Curve Name']=='CURVE 2'])
x122=np.array((x11_F5-float(x11_F5[:1])))
y122=np.array((y11_F5-float(y11_F5[:1])))


m2=(y122[-1])/x122[-1]
angle2 = np.degrees(np.arctan(m2))

rules1 = [x122[-1]>0, y122[-1]>0, m2>0]
rules2 = [x122[-1]>0, y122[-1]<0, m2<0]
rules3 = [x122[-1]<0, y122[-1]<0, m2>0]
rules4 = [x122[-1]<0, y122[-1]>0, m2<0]


if  all(rules1):
    theta2= 360-angle2
elif  all(rules2): 
    theta2= abs(angle2)
elif  all(rules3): 
    theta2= 180-angle2
elif  all(rules4): 
    theta2= 180-angle2
print(theta2)

x2_F5=x122*(math.cos(math.radians(theta2)))-y122*(math.sin(math.radians(theta2)))
y2_F5=x122*(math.sin(math.radians(theta2)))+y122*(math.cos(math.radians(theta2)))

f=interp1d(x2_F5,y2_F5)

distance = np.cumsum(np.sqrt( np.ediff1d(x2_F5, to_begin=0)**2 + np.ediff1d(y2_F5, to_begin=0)**2 ))
distance = distance/distance[-1]

fx2, fy2 = interp1d( distance, x2_F5 ), interp1d( distance, y2_F5 )

alpha = np.linspace(0, 1, 15)
x_regular2_F5, y_regular2_F5 = fx2(alpha), fy2(alpha)


#to calculate sinuosity
d2=np.zeros(14)
for i in range(14):
    d2[i]=np.sqrt((x_regular2_F5[i+1]-x_regular2_F5[i])**2+(y_regular2_F5[i+1]-y_regular2_F5[i])**2)
print(d2)

length2=sum(d2)
distance2=np.sqrt((x_regular2_F5[14]-x_regular2_F5[0])**2+(y_regular2_F5[14]-y_regular2_F5[0])**2)
sinuosity2_F5=length2/distance2
#---------------------------------------------------------
x12_F5=df_5['X-Coordinate '] [df_5['Curve Name']=='CURVE 3']
y12_F5=(df_5['Y-Coordinate '] [df_5['Curve Name']=='CURVE 3'])
x123=np.array((x12_F5-float(x12_F5[:1])))
y123=np.array((y12_F5-float(y12_F5[:1])))


m3=(y123[-1])/x123[-1]
angle3 = np.degrees(np.arctan(m3))

rules1 = [x123[-1]>0, y123[-1]>0, m3>0]
rules2 = [x123[-1]>0, y123[-1]<0, m3<0]
rules3 = [x123[-1]<0, y123[-1]<0, m3>0]
rules4 = [x123[-1]<0, y123[-1]>0, m3<0]


if  all(rules1):
    theta3= 360-angle3
elif  all(rules2): 
    theta3= abs(angle3)
elif  all(rules3): 
    theta3= 180-angle3
elif  all(rules4): 
    theta3= 180-angle3
print(theta3)

x3_F5=x123*(math.cos(math.radians(theta3)))-y123*(math.sin(math.radians(theta3)))
y3_F5=x123*(math.sin(math.radians(theta3)))+y123*(math.cos(math.radians(theta3)))


f=interp1d(x3_F5,y3_F5)
distance = np.cumsum(np.sqrt( np.ediff1d(x3_F5, to_begin=0)**2 + np.ediff1d(y3_F5, to_begin=0)**2 ))
distance = distance/distance[-1]

fx3, fy3 = interp1d( distance, x3_F5 ), interp1d( distance, y3_F5 )

alpha = np.linspace(0, 1, 15)
x_regular3_F5, y_regular3_F5 = fx3(alpha), fy3(alpha)


#to calculate sinuosity
d3=np.zeros(14)
for i in range(14):
    d3[i]=np.sqrt((x_regular3_F5[i+1]-x_regular3_F5[i])**2+(y_regular3_F5[i+1]-y_regular3_F5[i])**2)
print(d3)

length3=sum(d3)
distance3=np.sqrt((x_regular3_F5[14]-x_regular3_F5[0])**2+(y_regular3_F5[14]-y_regular3_F5[0])**2)
sinuosity3_F5=length3/distance3
#--------------------------------------------------
x13_F5=df_5['X-Coordinate '] [df_5['Curve Name']=='CURVE 4']
y13_F5=(df_5['Y-Coordinate '] [df_5['Curve Name']=='CURVE 4'])
x124=np.array((x13_F5-float(x13_F5[:1])))
y124=np.array((y13_F5-float(y13_F5[:1])))


m4=(y124[-1])/x124[-1]
angle4 = np.degrees(np.arctan(m4))

rules1 = [x124[-1]>0, y124[-1]>0, m4>0]
rules2 = [x124[-1]>0, y124[-1]<0, m4<0]
rules3 = [x124[-1]<0, y124[-1]<0, m4>0]
rules4 = [x124[-1]<0, y124[-1]>0, m4<0]


if  all(rules1):
    theta4= 360-angle4
elif  all(rules2): 
    theta4= abs(angle4)
elif  all(rules3): 
    theta4= 180-angle4
elif  all(rules4): 
    theta4= 180-angle4
print(theta4)


x4_F5=x124*(math.cos(math.radians(theta4)))-y124*(math.sin(math.radians(theta4)))
y4_F5=x124*(math.sin(math.radians(theta4)))+y124*(math.cos(math.radians(theta4)))


f=interp1d(x4_F5,y4_F5)
distance = np.cumsum(np.sqrt( np.ediff1d(x4_F5, to_begin=0)**2 + np.ediff1d(y4_F5, to_begin=0)**2 ))
distance = distance/distance[-1]

fx4, fy4 = interp1d( distance, x4_F5 ), interp1d( distance, y4_F5 )

alpha = np.linspace(0, 1, 15)
x_regular4_F5, y_regular4_F5 = fx4(alpha), fy4(alpha)


#to calculate sinuosity
d4=np.zeros(15)
for i in range(14):
    d4[i]=np.sqrt((x_regular4_F5[i+1]-x_regular4_F5[i])**2+(y_regular4_F5[i+1]-y_regular4_F5[i])**2)
print(d4)

length4=sum(d4)
distance4=np.sqrt((x_regular4_F5[14]-x_regular4_F5[0])**2+(y_regular4_F5[14]-y_regular4_F5[0])**2)
sinuosity4_F5=length4/distance4
#--------------------------------------------------
x14_F5=df_5['X-Coordinate '] [df_5['Curve Name']=='CURVE 5']
y14_F5=(df_5['Y-Coordinate '] [df_5['Curve Name']=='CURVE 5'])
x125=np.array((x14_F5-float(x14_F5[:1])))
y125=np.array((y14_F5-float(y14_F5[:1])))


m5=(y125[-1])/x125[-1]
angle5= np.degrees(np.arctan(m5))

rules1 = [x125[-1]>0, y125[-1]>0, m5>0]
rules2 = [x125[-1]>0, y125[-1]<0, m5<0]
rules3 = [x125[-1]<0, y125[-1]<0, m5>0]
rules4 = [x125[-1]<0, y125[-1]>0, m5<0]


if  all(rules1):
    theta5= 360-angle5
elif  all(rules2): 
    theta5= abs(angle5)
elif  all(rules3): 
    theta5= 180-angle5
elif  all(rules4): 
    theta5= 180-angle5


x5_F5=x125*(math.cos(math.radians(theta5)))-y125*(math.sin(math.radians(theta5)))
y5_F5=x125*(math.sin(math.radians(theta5)))+y125*(math.cos(math.radians(theta5)))

f=interp1d(x5_F5,y5_F5)
distance = np.cumsum(np.sqrt( np.ediff1d(x5_F5, to_begin=0)**2 + np.ediff1d(y5_F5, to_begin=0)**2 ))
distance = distance/distance[-1]

fx5, fy5 = interp1d( distance, x5_F5), interp1d( distance, y5_F5 )

alpha = np.linspace(0, 1, 15)
x_regular5_F5, y_regular5_F5 = fx5(alpha), fy5(alpha)


#to calculate sinuosity
d5=np.zeros(14)
for i in range(14):
    d5[i]=np.sqrt((x_regular5_F5[i+1]-x_regular5_F5[i])**2+(y_regular5_F5[i+1]-y_regular5_F5[i])**2)
print(d5)

length5=sum(d5)
distance5=np.sqrt((x_regular5_F5[14]-x_regular5_F5[0])**2+(y_regular5_F5[14]-y_regular5_F5[0])**2)
sinuosity5_F5=length5/distance5

#center point##################################
x_all=np.zeros(5)
x_all[0]=x1_F5.iloc[0]
x_all[1]=x11_F5.iloc[0]
x_all[2]=x12_F5.iloc[0]
x_all[3]=x13_F5.iloc[0]
x_all[4]=x14_F5.iloc[0]
x_disc_F5=x_all.mean()
y_all=np.zeros(5)
y_all[0]=y1_F5.iloc[0]
y_all[1]=y11_F5.iloc[0]
y_all[2]=y12_F5.iloc[0]
y_all[3]=y13_F5.iloc[0]
y_all[4]=y14_F5.iloc[0]
y_disc_F5=y_all.mean()

###the last point of curves in F5
y1_tip_F5=y1_F5.iloc[-1]
x1_tip_F5=x1_F5.iloc[-1]
y2_tip_F5=y11_F5.iloc[-1]
x2_tip_F5=x11_F5.iloc[-1]
y3_tip_F5=y12_F5.iloc[-1]
x3_tip_F5=x12_F5.iloc[-1]
y4_tip_F5=y13_F5.iloc[-1]
x4_tip_F5=x13_F5.iloc[-1]
y5_tip_F5=y14_F5.iloc[-1]
x5_tip_F5=x14_F5.iloc[-1]

###the first point of curves in F5
y1_base_F5=y1_F4.iloc[0]
x1_base_F5=x1_F4.iloc[0]
y2_base_F5=y11_F4.iloc[0]
x2_base_F5=x11_F4.iloc[0]
y3_base_F5=y12_F4.iloc[0]
x3_base_F5=x12_F4.iloc[0]
y4_base_F5=y13_F4.iloc[0]
x4_base_F5=x13_F4.iloc[0]
y5_base_F5=y14_F4.iloc[0]
x5_base_F5=x14_F4.iloc[0]



###################################################################### Arm angle F5

##disc new= disc points-disc points

x_disc_F5_n=x_disc_F5-x_disc_F5
y_disc_F5_n=y_disc_F5-y_disc_F5



## basal point new= basal point - disc point

y1_base_F5_n=y1_F5.iloc[0]-y_disc_F5
x1_base_F5_n=x1_F5.iloc[0]-x_disc_F5
y2_base_F5_n=y11_F5.iloc[0]-y_disc_F5
x2_base_F5_n=x11_F5.iloc[0]-x_disc_F5
y3_base_F5_n=y12_F5.iloc[0]-y_disc_F5
x3_base_F5_n=x12_F5.iloc[0]-x_disc_F5
y4_base_F5_n=y13_F5.iloc[0]-y_disc_F5
x4_base_F5_n=x13_F5.iloc[0]-x_disc_F5
y5_base_F5_n=y14_F5.iloc[0]-y_disc_F5
x5_base_F5_n=x14_F5.iloc[0]-x_disc_F5

x_base=[x1_base_F5_n, x2_base_F5_n, x3_base_F5_n, x4_base_F5_n, x5_base_F5_n]
y_base=[y1_base_F5_n, y2_base_F5_n, y3_base_F5_n, y4_base_F5_n, y5_base_F5_n]

##tip point new= tip point - disc point

y1_tip_F5_n=y1_F5.iloc[-1]-y_disc_F5
x1_tip_F5_n=x1_F5.iloc[-1]-x_disc_F5
y2_tip_F5_n=y11_F5.iloc[-1]-y_disc_F5
x2_tip_F5_n=x11_F5.iloc[-1]-x_disc_F5
y3_tip_F5_n=y12_F5.iloc[-1]-y_disc_F5
x3_tip_F5_n=x12_F5.iloc[-1]-x_disc_F5
y4_tip_F5_n=y13_F5.iloc[-1]-y_disc_F5
x4_tip_F5_n=x13_F5.iloc[-1]-x_disc_F5
y5_tip_F5_n=y14_F5.iloc[-1]-y_disc_F5
x5_tip_F5_n=x14_F5.iloc[-1]-x_disc_F5

x_tip=[x1_tip_F5_n, x2_tip_F5_n, x3_tip_F5_n, x4_tip_F5_n, x5_tip_F5_n]
y_tip=[y1_tip_F5_n, y2_tip_F5_n, y3_tip_F5_n, y4_tip_F5_n, y5_tip_F5_n]

## alpha is the angle of disc center to basal point

m_alpha1= y1_base_F5_n/x1_base_F5_n
m_alpha2= y2_base_F5_n/x2_base_F5_n
m_alpha3= y3_base_F5_n/x3_base_F5_n
m_alpha4= y4_base_F5_n/x4_base_F5_n
m_alpha5= y5_base_F5_n/x5_base_F5_n

alpha1= np.degrees(np.arctan(m_alpha1))
alpha2= np.degrees(np.arctan(m_alpha2))
alpha3= np.degrees(np.arctan(m_alpha3))
alpha4= np.degrees(np.arctan(m_alpha4))
alpha5= np.degrees(np.arctan(m_alpha5))

alpha_total=[alpha1, alpha2, alpha3, alpha4, alpha5]


## betha is the angle of disc center to tip point

m_betha1= y1_tip_F5_n/x1_tip_F5_n
m_betha2= y2_tip_F5_n/x2_tip_F5_n
m_betha3= y3_tip_F5_n/x3_tip_F5_n
m_betha4= y4_tip_F5_n/x4_tip_F5_n
m_betha5= y5_tip_F5_n/x5_tip_F5_n

betha1= np.degrees(np.arctan(m_betha1))
betha2= np.degrees(np.arctan(m_betha2))
betha3= np.degrees(np.arctan(m_betha3))
betha4= np.degrees(np.arctan(m_betha4))
betha5= np.degrees(np.arctan(m_betha5))

betha_total= [betha1, betha2, betha3, betha4, betha5]
#########################################################

alpha_total_n=np.zeros(5)

for j in range (5):
    if x_base[j]<0 and y_base[j]>0:
        alpha_total_n[j]= 180+alpha_total[j]
    
    elif x_base[j]<0 and y_base[j]<0:
        alpha_total_n[j]= 180-alpha_total[j]
    
    elif  x_base[j]>0 and y_base[j]<0:
        alpha_total_n[j]= -alpha_total[j]
        
    elif x_base[j]>0 and y_base[j]>0:
        alpha_total_n[j]=alpha_total[j]

print (alpha_total_n)
#######################################################
betha_total_n=np.zeros(5)
   
for k in range (5):
    if x_tip[k]<0 and y_tip[k]>0:
        betha_total_n[k]= 180+betha_total[k]
    
    elif x_tip[k]<0 and y_tip[k]<0:
        betha_total_n[k]= 180-betha_total[k]
    
    elif x_tip[k]>0 and y_tip[k]<0:
        betha_total_n[k]= -betha_total[k]
        
    elif x_tip[k]>0 and y_tip[k]>0:
        betha_total_n[k]=betha_total[k]
        
print(betha_total_n)

######################################################

#### the difference between two lines

phi1= betha_total_n[0]-alpha_total_n[0]
phi2= betha_total_n[1]-alpha_total_n[1]
phi3= betha_total_n[2]-alpha_total_n[2]
phi4= betha_total_n[3]-alpha_total_n[3]
phi5= betha_total_n[4]-alpha_total_n[4]

phi_total_F5= [phi1, phi2, phi3, phi4, phi5]
phi_total_F5

##################################################

for i in range(5):
    if phi_total_F5[i] > 0:
        if phi_total_F5[i] < 180:
            print ('phi positive') 
        else:
            print('phi negative') 
    else:
        if abs(phi_total_F5[i]) <180:
            print('phi negative') 
        else:
            print('phi positive')    
            
##### ASSUMPTION: ABS PHI <180


# ## Variable 3: Disc displacement of disc center between time frames 1-5

# In[ ]:


##distance changes of disc center time frames1-2 , 2-3, 3-4, 4-5
disc_dist_F1_F2=np.sqrt((x_disc_F2-x_disc_F1)**2+(y_disc_F2-y_disc_F1)**2)
disc_dist_F2_F3=np.sqrt((x_disc_F3-x_disc_F2)**2+(y_disc_F3-y_disc_F2)**2)
disc_dist_F3_F4=np.sqrt((x_disc_F4-x_disc_F3)**2+(y_disc_F4-y_disc_F3)**2)
disc_dist_F4_F5=np.sqrt((x_disc_F5-x_disc_F4)**2+(abs(y_disc_F5)-y_disc_F4)**2)


# ## Variable 4: Arm slip angle between time frames 1-5

# In[ ]:



################################################################to calculate the anglies between two lines frame 1-2
y_disc_F1_n=y_disc_F1-y_disc_F1
x_disc_F1_n=x_disc_F1-x_disc_F1
y1_base_F1_n=y1_base_F1-y_disc_F1
x1_base_F1_n=x1_base_F1-x_disc_F1
y_disc_F2_n=y_disc_F2-y_disc_F1
x_disc_F2_n=x_disc_F2-x_disc_F1

m1=(y1_base_F1_n-y_disc_F1_n)/(x1_base_F1_n-x_disc_F1_n) #=m_disc_base_arm1
m2=(y_disc_F2_n-y_disc_F1_n)/(x_disc_F2_n-x_disc_F1_n)
g=abs((m2-m1)/(1+(m1*m2)))
angle_ = np.degrees(np.arctan(g))

y2_base_F1_n=y2_base_F1-y_disc_F1
x2_base_F1_n=x2_base_F1-x_disc_F1

m11=(y2_base_F1_n-x_disc_F1_n)/(x2_base_F1_n-x_disc_F1_n)   ##=m_disc_base_arm2
m21=(y_disc_F2_n-x_disc_F1_n)/(x_disc_F2_n-x_disc_F1_n)
g2=abs((m21-m11)/(1+(m11*m21)))
g2
if x_disc_F2_n<0:
    angle_2 = 180-np.degrees(np.arctan(g2))
else:
    angle_2=np.degrees(np.arctan(g2))

y3_base_F1_n=y3_base_F1-y_disc_F1
x3_base_F1_n=x3_base_F1-x_disc_F1

m12=(y3_base_F1_n-y_disc_F1_n)/(x3_base_F1_n-x_disc_F1_n) #=m_disc_base_arm3
m22=(y_disc_F2_n-y_disc_F1_n)/(x_disc_F2_n-x_disc_F1_n)
g3=abs((m22-m12)/(1+(m12*m22)))
g3
angle_3 = 180-np.degrees(np.arctan(g3))

y4_base_F1_n=y4_base_F1-y_disc_F1
x4_base_F1_n=x4_base_F1-x_disc_F1

m13=(y4_base_F1_n-y_disc_F1_n)/(x4_base_F1_n-x_disc_F1_n) #=m_disc_base_arm4
m23=(y_disc_F2_n-y_disc_F1_n)/(x_disc_F2_n-x_disc_F1_n)
g4=abs((m23-m13)/(1+(m13*m23)))
g4
angle_4 = 180-np.degrees(np.arctan(g4))

y5_base_F1_n=y5_base_F1-y_disc_F1
x5_base_F1_n=x5_base_F1-x_disc_F1

m14=(y5_base_F1_n-y_disc_F1_n)/(x5_base_F1_n-x_disc_F1_n) #=m_disc_base_arm5
m24=(y_disc_F2_n-y_disc_F1_n)/(x_disc_F2_n-x_disc_F1_n)
g5=abs((m24-m14)/(1+(m14*m24)))
g5
if x_disc_F2_n>0:
    angle_5=180-np.degrees(np.arctan(g5))
else:
    angle_5 = np.degrees(np.arctan(g5))
Angles_F1_F2=[angle_, angle_2,angle_3,angle_4, angle_5]


############################################################to calculate the anglies between two lines frame 2-3
y_disc_F2_n=y_disc_F2-y_disc_F2
x_disc_F2_n=x_disc_F2-x_disc_F2
y1_base_F2_n=y1_base_F2-y_disc_F2
x1_base_F2_n=x1_base_F2-x_disc_F2
y_disc_F3_n=y_disc_F3-y_disc_F2
x_disc_F3_n=x_disc_F3-x_disc_F2


m1=(y1_base_F2_n-y_disc_F2_n)/(x1_base_F2_n-x_disc_F2_n) 
m2=(y_disc_F3_n-y_disc_F2_n)/(x_disc_F3_n-x_disc_F2_n)
g=abs((m2-m1)/(1+(m1*m2)))
angle_ = np.degrees(np.arctan(g))

y2_base_F2_n=y2_base_F2-y_disc_F2
x2_base_F2_n=x2_base_F2-x_disc_F2

m11=(y2_base_F2_n-y_disc_F2_n)/(x2_base_F2_n-x_disc_F2_n) 
m21=(y_disc_F3_n-y_disc_F2_n)/(x_disc_F3_n-x_disc_F2_n)
g2=abs((m21-m11)/(1+(m11*m21)))
g2

if x_disc_F3_n<0:
    angle_2=180-np.degrees(np.arctan(g2))
else:
    angle_2 = np.degrees(np.arctan(g2))

y3_base_F2_n=y3_base_F2-y_disc_F2
x3_base_F2_n=x3_base_F2-x_disc_F2

m12=(y3_base_F2_n-y_disc_F2_n)/(x3_base_F2_n-x_disc_F2_n) 
m22=(y_disc_F3_n-y_disc_F2_n)/(x_disc_F3_n-x_disc_F2_n)
g3=abs((m22-m12)/(1+(m12*m22)))
g3
angle_3 = 180-np.degrees(np.arctan(g3))

y4_base_F2_n=y4_base_F2-y_disc_F2
x4_base_F2_n=x4_base_F2-x_disc_F2

m13=(y4_base_F2_n-y_disc_F2_n)/(x4_base_F2_n-x_disc_F2_n) 
m23=(y_disc_F3_n-y_disc_F2_n)/(x_disc_F3_n-x_disc_F2_n)
g4=abs((m23-m13)/(1+(m13*m23)))
g4
angle_4 = 180-np.degrees(np.arctan(g4))

y5_base_F2_n=y5_base_F2-y_disc_F2
x5_base_F2_n=x5_base_F2-x_disc_F2

m14=(y5_base_F2_n-y_disc_F2_n)/(x5_base_F2_n-x_disc_F2_n) 
m24=(y_disc_F3_n-y_disc_F2_n)/(x_disc_F3_n-x_disc_F2_n)
g5=abs((m24-m14)/(1+(m14*m24)))
g5
if x_disc_F3_n>0:
    angle_5=180-np.degrees(np.arctan(g5))
else:
    angle_5 = np.degrees(np.arctan(g5))
Angles_F2_F3=[angle_, angle_2,angle_3,angle_4, angle_5]

#################################################################to calculate the anglies between two lines frame 3-4
y_disc_F3_n=(y_disc_F3-y_disc_F3)
x_disc_F3_n=(x_disc_F3-x_disc_F3)
y1_base_F3_n=(y1_base_F3-y_disc_F3)
x1_base_F3_n=(x1_base_F3-x_disc_F3)
y_disc_F4_n=(y_disc_F4-y_disc_F3)
x_disc_F4_n=(x_disc_F4-x_disc_F3)


m1=(y1_base_F3_n-y_disc_F3_n)/(x1_base_F3_n-x_disc_F3_n) 
m2=(y_disc_F4_n-y_disc_F3_n)/(x_disc_F4_n-x_disc_F3_n)
g=abs((m2-m1)/(1+(m1*m2)))
angle_ = np.degrees(np.arctan(g))

y2_base_F3_n=y2_base_F3-y_disc_F3
x2_base_F3_n=x2_base_F3-x_disc_F3

m11=(y2_base_F3_n-y_disc_F3_n)/(x2_base_F3_n-x_disc_F3_n) 
m21=(y_disc_F4_n-y_disc_F3_n)/(x_disc_F4_n-x_disc_F3_n)
g2=abs((m21-m11)/(1+(m11*m21)))
g2

if x_disc_F4_n<0:
    angle_2=180-np.degrees(np.arctan(g2))
else:
    angle_2 = np.degrees(np.arctan(g2))

y3_base_F3_n=y3_base_F3-y_disc_F3
x3_base_F3_n=x3_base_F3-x_disc_F3

m12=(y3_base_F3_n-y_disc_F3_n)/(x3_base_F3_n-x_disc_F3_n) 
m22=(y_disc_F4_n-y_disc_F3_n)/(x_disc_F4_n-x_disc_F3_n)
g3=abs((m22-m12)/(1+(m12*m22)))
g3
angle_3 = 180-np.degrees(np.arctan(g3))

y4_base_F3_n=y4_base_F3-y_disc_F3
x4_base_F3_n=x4_base_F3-x_disc_F3

m13=(y4_base_F3_n-y_disc_F3_n)/(x4_base_F3_n-x_disc_F3_n) 
m23=(y_disc_F4_n-y_disc_F3_n)/(x_disc_F4_n-x_disc_F3_n)
g4=abs((m23-m13)/(1+(m13*m23)))
g4
angle_4 = 180-np.degrees(np.arctan(g4))

y5_base_F3_n=y5_base_F3-y_disc_F3
x5_base_F3_n=x5_base_F3-x_disc_F3

m14=(y5_base_F3_n-y_disc_F3_n)/float(x5_base_F3_n-x_disc_F3_n) 
m24=(y_disc_F4_n-y_disc_F3_n)/(x_disc_F4_n-x_disc_F3_n)
g5=abs((m24-m14)/(1+(m14*m24)))
g5
if x_disc_F4_n>0:
    angle_5=180-np.degrees(np.arctan(g5))
else:
    angle_5 = np.degrees(np.arctan(g5))
Angles_F3_F4=[angle_, angle_2,angle_3,angle_4, angle_5]

############################################################### to calculate the anglies between two lines frame 4-5
y_disc_F4_n=y_disc_F4-y_disc_F4
x_disc_F4_n=x_disc_F4-x_disc_F4
y1_base_F4_n=y1_base_F4-y_disc_F4
x1_base_F4_n=x1_base_F4-x_disc_F4
y_disc_F5_n=y_disc_F5-y_disc_F4
x_disc_F5_n=x_disc_F5-x_disc_F4


m1=(y1_base_F4_n-y_disc_F4_n)/(x1_base_F4_n-x_disc_F4_n) 
m2=(y_disc_F5_n-y_disc_F4_n)/(x_disc_F5_n-x_disc_F4_n)
g=abs((m2-m1)/(1+(m1*m2)))
angle_ = np.degrees(np.arctan(g))

y2_base_F4_n=y2_base_F4-y_disc_F4
x2_base_F4_n=x2_base_F4-x_disc_F4

m11=(y2_base_F4_n-y_disc_F4_n)/(x2_base_F4_n-x_disc_F4_n) 
m21=(y_disc_F5_n-y_disc_F4_n)/(x_disc_F5_n-x_disc_F4_n)
g2=abs((m21-m11)/(1+(m11*m21)))
g2

if x_disc_F5_n<0:
    angle_5=180-np.degrees(np.arctan(g2))
else:
    angle_5 = np.degrees(np.arctan(g2))

y3_base_F4_n=y3_base_F4-y_disc_F4
x3_base_F4_n=x3_base_F4-x_disc_F4

m12=(y3_base_F4_n-y_disc_F4_n)/(x3_base_F4_n-x_disc_F4_n) 
m22=(y_disc_F5_n-y_disc_F4_n)/(x_disc_F5_n-x_disc_F4_n)
g3=abs((m22-m12)/(1+(m12*m22)))
g3
angle_3 = 180-np.degrees(np.arctan(g3))

y4_base_F4_n=y4_base_F4-y_disc_F4
x4_base_F4_n=x4_base_F4-x_disc_F4

m13=(y4_base_F4_n-y_disc_F4_n)/(x4_base_F4_n-x_disc_F4_n) 
m23=(y_disc_F5_n-y_disc_F4_n)/(x_disc_F5_n-x_disc_F4_n)
g4=abs((m23-m13)/(1+(m13*m23)))
g4
angle_4 = 180-np.degrees(np.arctan(g4))

y5_base_F4_n=y5_base_F4-y_disc_F4
x5_base_F4_n=x5_base_F4-x_disc_F4

m14=(y5_base_F4_n-y_disc_F4_n)/(x5_base_F4_n-x_disc_F4_n) 
m24=(y_disc_F5_n-y_disc_F4_n)/(x_disc_F5_n-x_disc_F4_n)
g5 = abs((m24-m14)/(1+(m14*m24)))

if x_disc_F5_n>0:
    angle_5=180-np.degrees(np.arctan(g5))
else:
    angle_5 = np.degrees(np.arctan(g5))
Angles_F4_F5=[angle_, angle_2,angle_3,angle_4, angle_5]


# Mona.goharimanesh@ugent.be/ July, 2022
