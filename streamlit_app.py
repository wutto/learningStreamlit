###################################################################
##[0] Loading required modules (standard modules)
import math, os, copy

##Modules for Data science
import pandas as pd
import numpy as np
import pickle

##Module for Visualization
import streamlit as st
import streamlit.components.v1 as components
#import plotly.express as px
#import plotly.graph_objects as go

##Module for machine learning
#from sklearn.linear_model import LinearRegression
###################################################################
colorSet1 = {'NCBI': '#5b86ae', 'JGI': '#85bc37'}
colorSet2 = {"NRPS":"#2D4D42","PKS":"#718878","RiPP":"#C4D0CC","other":"#DED6D3","terpene":"#B6AEAB"}
colorSet2_1 = {"NRPS":"#95A1AE","PKS":"#C8CFD6","RiPP":"#919F89","other":"#EFBC68","terpene":"#FFB7A1"}
colorSet3 = {"NRPS-like":"#4D4138","NRPS":"#6D5C4F","CDPS":"#917A69","NAPAA":"#B5A599","thioamide-NRP":"#D5CCC5","siderophore":"#21361E","betalactone":"#2C4828","indole":"#304E2C","other":"#3B6036","phosphonate":"#45713F","arylpolyene":"#52854B","hserlactone":"#5F9B57","butyrolactone":"#71AB69","resorcinol":"#8CBC86","ectoine":"#ABCEA6","ladderane":"#D9E9D7","T3PKS":"#6A210A","T1PKS":"#BA3A12","transAT-PKS":"#E94C1B","transAT-PKS-like":"#F08766","prodigiosin":"#F4A890","PKS-like":"#FAD5CA","borosin_MTYorD":"#091519","fungal-RiPP":"#132C35","microviridin":"#1E4452","RiPP-like":"#285C6E","cyclic-lactone-autoinducer":"#347890","lanthipeptide-class-iv":"#3E8FAC","RRE-containing":"#4DA1BF","thiopeptide":"#8EC3D6","linaridin":"#B3D6E3","lanthipeptide-class-ii":"#D5E8EF","terpene":"#FABE00"
}
##Define function
def genSankey(df,cat_cols=[],value_cols='',title='Sankey Diagram'):
    # maximum of 6 value cols -> 6 colors
    colorPalette = ['#9DAAA2','#8DD1C5','#EDCC8B','#E8B298','#D4A29C','#A26360']##color from https://www.pinterest.ch/pin/39688040457483199/
    labelList = []
    colorNumList = []
    for catCol in cat_cols:
        labelListTemp =  list(set(df[catCol].values))
        colorNumList.append(len(labelListTemp))
        labelList = labelList + labelListTemp
        
    # remove duplicates from labelList
    labelList = list(dict.fromkeys(labelList))
    
    # define colors based on number of levels
    colorList = []
    for idx, colorNum in enumerate(colorNumList):
        colorList = colorList + [colorPalette[idx]]*colorNum
        
    # transform df into a source-target pair
    for i in range(len(cat_cols)-1):
        if i==0:
            sourceTargetDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            sourceTargetDf.columns = ['source','target','count']
        else:
            tempDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            tempDf.columns = ['source','target','count']
            sourceTargetDf = pd.concat([sourceTargetDf,tempDf])
        sourceTargetDf = sourceTargetDf.groupby(['source','target']).agg({'count':'sum'}).reset_index()
        
    # add index for source-target pair
    sourceTargetDf['sourceID'] = sourceTargetDf['source'].apply(lambda x: labelList.index(x))
    sourceTargetDf['targetID'] = sourceTargetDf['target'].apply(lambda x: labelList.index(x))
    
    # creating the sankey diagram
    data = dict(
        type='sankey',
        node = dict(
          pad = 6,
          thickness = 50,
          line = dict(
            color = "black",
            width = 0.9
          ),
          label = labelList,
          color = colorList
        ),
        link = dict(source = sourceTargetDf['sourceID'],target = sourceTargetDf['targetID'],value = sourceTargetDf['count']))
    layout =  dict(title = title,font = dict(size = 60))
    fig = dict(data=[data], layout=layout)
    return fig
#############################################
def createBarPercentCategory(df_bgc,df_genome,rank):
    df_merge = pd.merge(df_bgc,df_genome, how="left", on=["path"])
    df_plot = df_merge.pivot_table(values='count', index = ["category",rank], aggfunc=np.sum).reset_index()
    data = df_plot.values.tolist()
    taxoCount = {}
    for record in data:
        key = record[1]
        if key not in taxoCount:
            value = 0
        else:
            value = taxoCount[key]
        value = value + record[2]##concatenate count
        taxoCount[key] = value##update value in dictionary
    percent = []
    for record in data:
        record.append(record[-1]/taxoCount[record[1]]*100.0)
        percent.append(record)

    df_plot_percent = pd.DataFrame(percent, columns=["category",rank,"count","percent"])
    return df_plot_percent
#############################################
def createBarPercentProduct(df_bgc,df_genome,rank):
    df_merge = pd.merge(df_bgc,df_genome, how="left", on=["path"])
    df_plot = df_merge.pivot_table(values='count', index = ["product",rank], aggfunc=np.sum).reset_index()
    data = df_plot.values.tolist()
    taxoCount = {}
    for record in data:
        key = record[1]
        if key not in taxoCount:
            value = 0
        else:
            value = taxoCount[key]
        value = value + record[2]##multiply count
        taxoCount[key] = value##update value in dictionary
    percent = []
    for record in data:
        record.append(record[-1]/taxoCount[record[1]]*100.0)
        percent.append(record)

    df_plot_percent = pd.DataFrame(percent, columns=["product",rank,"count","percent"])
    return df_plot_percent
#############################################
def createRadarCategory(df_bgc,df_genome,rank):
    df_merge = pd.merge(df_bgc,df_genome, how="left", on=["path"])
    df_plot = df_merge.pivot_table(values='count', index = ["category",rank], aggfunc=np.sum).reset_index()
    data = df_plot.values.tolist()
    
    ##Calculate percentage and average BGC per genome
    taxoCount = {}
    for record in data:
        key = record[1]
        if key not in taxoCount:
            value = 0
        else:
            value = taxoCount[key]
        value = value + record[2]##multiply count
        taxoCount[key] = value##update value in dictionary
        
    taxoPerGenome = {}
    for record in df_genome.filter(items=['path', rank]).values.tolist():
        key = record[1]
        if key not in taxoPerGenome:
            value = 0
        else:
            value = taxoPerGenome[key]
        value = value + 1 ##count number of genome in each taxon
        taxoPerGenome[key] = value##update value in dictionary
        
    percent = []
    for record in data:
        record.append(record[2]/taxoCount[record[1]]*100.0)
        record.append(record[2]/taxoPerGenome[record[1]])
        #record.append(taxoPerGenome[record[1]])
        percent.append(record)
        
    df_plot_percent = pd.DataFrame(percent, columns=["category",rank,"count","percent","pergenome"])
    return df_plot_percent
#############################################
def createRadarProduct(df_bgc,df_genome,rank,memberList):
    df_merge = pd.merge(df_bgc,df_genome, how="left", on=["path"])
    df_plot = df_merge.pivot_table(values='count', index = ["category","product",rank], aggfunc=np.sum).reset_index()
    df_plot = df_plot[df_plot[rank].isin(memberList)]
    data = df_plot.values.tolist()
    
    ##Calculate percentage and average BGC per genome
    taxoCount = {}
    for record in data:
        key = record[2]
        if key not in taxoCount:
            value = 0
        else:
            value = taxoCount[key]
        value = value + record[3]##multiply count
        taxoCount[key] = value##update value in dictionary
        
    taxoPerGenome = {}
    for record in df_genome.filter(items=['path', rank]).values.tolist():
        key = record[1]
        if key not in taxoPerGenome:
            value = 0
        else:
            value = taxoPerGenome[key]
        value = value + 1 ##count number of genome in each taxon
        taxoPerGenome[key] = value##update value in dictionary
        
    percent = []
    for record in data:
        record.append(record[3]/taxoCount[record[2]]*100.0)
        record.append(record[3]/taxoPerGenome[record[2]])
        percent.append(record)
        
    df_plot_percent = pd.DataFrame(percent, columns=["category","product",rank,"count","percent","pergenome"])
    return df_plot_percent

def filterDf(df,season,temp,tempRange,humidity,humidityRange):
    tempMin = temp - int(tempRange.split(" ")[-1])
    tempMax = temp + int(tempRange.split(" ")[-1])
    humidityMin = humidity - int(humidityRange.split(" ")[-1])
    humidityMax = humidity + int(humidityRange.split(" ")[-1])
    filtered_df = df[(df["season"] == season) & (df["optTemp"] >= tempMin) & (df["optTemp"] <= tempMax) & (df["humidity"] >= humidityMin) & (df["humidity"] <= humidityMax)]
    return filtered_df
###################################################################
##[1] Loading input file
##Summary and visualize GBK feature
tableEnv = "./environment.xlsx"
dfEnv = pd.read_excel(tableEnv, sheet_name="summary")
###################################################################
##[2] Visualizing data
##Set layout
st.set_page_config(layout="wide")
st.header("Pathogen prediction")
st.write("-----"*100)

#Preparing data frame
#st.dataframe(df_genome, use_container_width=True)
#st.write("-----"*100)

#df_bgc_sum = df_bgc.pivot_table(values='count', index = ["path"], aggfunc=np.sum).reset_index()
#df_genome_bgc = pd.merge(df_bgc_sum,df_genome, how="left", on=["path"])

#st.dataframe(df_genome_bgc, use_container_width=True)
#st.write("-----"*100)
####################################################
##Options
    ##Season
season = st.selectbox("Please select the season ",("rainy","winter"), index = 0)
    
#Set optimum temp
col_temp, col_tempRange = st.columns([8,2])
with col_temp:
    temp = st.slider('Please select the optimum temperature (Celsius)', 10, 50, 32)
with col_tempRange:
    tempRange = st.selectbox("Range",(" 0","+/- 1","+/- 5","+/- 10"), index = 2)

#Set humidity
humidity_col, humidityRange_col = st.columns([8,2])
with humidity_col:
    humidity = st.slider('Please select the humidity (%)', 10, 100, 80)
with humidityRange_col:
    humidityRange = st.selectbox("Range ",(" 0","+/- 1","+/- 5","+/- 10"), index = 2)
    
#Displaying summary of input
#st.write("-----"*100)
#st.write("temp:",temp)
#st.write("tempRange:",tempRange)
#st.write("humidity:",humidity)
#st.write("humidityRange:",humidityRange)
filterDf = filterDf(dfEnv,season,temp,tempRange,humidity,humidityRange)

if len(filterDf) > 0:
    ##Plot barchart | abundance of pathogen per taxonomy
    filterDf["combined"] = filterDf["kingdom"] + "_" + filterDf["compartment"]
    fig = px.bar(filterDf, x="species", y="abundance",
             color="combined",
             color_discrete_map = colorSet1,
             height=600,
             width=1800,
             labels={'abundance':'Percent abundance of species'},
             title="Summary number of pathogen")

    fig.update_traces(marker_line_color = 'black', marker_line_width = 1, opacity = 1)
    st.plotly_chart(fig, use_container_width=True)
    st.write("-----"*50)

    ##Plot bubble
    figBubble = px.scatter(filterDf, x="optTemp", y="humidity",size="abundance", color="kingdom", symbol="compartment", height=800, hover_name="species", size_max=50)
    figBubble.update_layout(
    title='Bubble plot: Optemp vs. Humidity',
    xaxis=dict(title='Optimum temperature',gridcolor='white',gridwidth=2),
    yaxis=dict(title='Humidity',gridcolor='white',gridwidth=2)
    )

    st.plotly_chart(figBubble, use_container_width=True)
    st.write("** Size of bubbles represent percent abundance")
    st.write("-----"*100)


    ##Plot sankey diagram
    dataSankey = genSankey(filterDf,cat_cols=["compartment","kingdom","phylum","class","family","species"],value_cols='abundance',title='Percent abundance')
    figSankey = go.Figure(dataSankey)
    figSankey.update_layout(title_text="Percent abundance", font_size=12)
    figSankey.update_layout(height=1000, width=1200)
    st.plotly_chart(figSankey, use_container_width=True)


##End
