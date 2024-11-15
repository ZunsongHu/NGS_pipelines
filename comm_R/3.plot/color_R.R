# Colors -------------------
cols_in_default=c('#3283FE','#F0A0FF','#0075DC','#5A5156','#E4E1E3','#F6222E','#FE00FA','#16FF32','#3283FE','#FEAF16','#B00068','#1CFFCE','#90AD1C','#2ED9FF',
                  '#DEA0FD','#AA0DFE','#F8A19F','#325A9B','#C4451C','#1C8356','#85660D','#B10DA1','#FBE426','#1CBE4F','#FA0087',
                  '#FC1CBF','#F7E1A0','#C075A6','#782AB6','#AAF400','#BDCDFF','#822E1C','#B5EFB5','#7ED7D1','#1C7F93','#D85FF7',
                  '#683B79','#66B0FF','#3B00FB')

cols_pool=c('red4','#FFA405','#2ED9FF','#B10DA1','#85660D','#325A9B','#1CBE4F','lightslateblue',
            '#3283FE','#7ED7D0','#683B79','#F8A19F','#90AD1C','#1C8356','#F7E1A0',
            '#822E1C','#B5EFB5','#BDCDFF','#1C7F93','#FBE426','#1CFFCE','#16FF32')

"OtherVDJ"='#3B00FB'
cols_in_default=cols_pool
cols_pool1=c('#1F78B5','pink','red4','skyblue','seagreen2','orangered','#00FF00','#1E90FF','#3E9F32','#469990','#66C2A6',
             '#808000','#A8DD00','black','blue3','#CCCC33','lightslateblue',
             'cyan','#D27B1C86','darkgoldenrod4','#E6BEFF','#FF00B6','#FFA620','gold2','grey40','grey75',
             'magenta3')

cols_2=c('#CA6D67','#68B6BA')
cols_4=c('red4','skyblue','#CCCC33','magenta3')
cols_5=c('#DE6E67','#717335',"#2AB47A",'#34ACE2',"#AF70A5")
cols_15=c('red4','#FFA405','#2ED9FF','#B10DA1','#85660D','#325A9B','#1CBE4F','lightslateblue',
          '#3283FE','#7ED7D0','#683B79','#F8A19F','#90AD1C','#1C8356','#F7E1A0')


cols_66=c('#FF33CC','#6633CC','#0000FF','#800000','#FF6633','#FFFF00','#FF0000','#008000','#FF00FF','#00FFFF','#FFA500','#EE82EE','#00FF00','#4B0082','#008080',
          '#808000','#800080','#008080','#C0C0C0','#808080','#999999','#333333','#666666','#CCCCCC','#663399','#FF99CC','#FFCC99','#FFFF99','#CCFF99','#000000',
          '#99FFCC','#99CCFF','#CC99FF','#FF99FF','#66CCCC','#99CCCC','#CCCC99','#CC99CC','#FF6699','#FF9966','#FFFF66','#99FF66','#66FF99','#66CCFF','#A52A2A',
          '#9966CC','#FF66CC','#339966','#66CC66','#CCCC66','#CC66CC','#FF3333','#FFFF33','#99FF33','#33FF66','#33CCFF',
          '#006633','#336666','#666633','#663366','#FF0033','#FF3300','#FFFF00','#99FF00','#00FF33','#0099FF')





c('#F0A0FF','#0075DC','#993F00','#4C005C','#005C31','#2BCE48','#191919',
  '#FFCC99','#808080','#94FFB5','#8F7C00','#9DCC00','#C20088','#003380',
  '#FFA8BB','#426600','#FF0010','#5EF1F2','#00998F')




c('#B00068','#FA0087','#FC1CBF','#C075A6','#782AB6','#AAF400',
  '#D85FF7','#FE00FA','#FEAF16','#66B0FF',
  '#DEA0FD','#AA0DFE','#C4451C','#5A5156','#E4E1E3',
  
  '#3B00FB')


col_default=c('#F6222E','#FFA405','#1C8356','#B10DA1','#85660D','#1CFFCE','#FBE426',
              '#F8A19F','#683B79','#3283FE','#16FF32','#7ED7D0','#90AD1C','#822E1C',
              '#F7E1A0','#2ED9FF','#B5EFB5','#BDCDFF','#1C7F93','#1CBE4F','#325A9B')



col_panel=c('#AAFFC3','#800000','#9A6324','#F032E6','#E6194B','#F58231','#4363D8',
            'gold2','pink','#1F78B5','magenta3','#A8DD00',
            '#66C2A6','seagreen2','black','skyblue','orangered','orangered','#3E9F32',
            '#1E90FF','blue3','red4','darkgoldenrod4',
            '#FFA620','red4','lightslateblue','#CCCC33','cyan'
)

# color -------------------
#BALL color




subtypeCol=c() 
{
  subtypeCol["ETV6::RUNX1"]="gold2"
  subtypeCol["ETV6::RUNX1-like"]="pink"
  subtypeCol["ETV6::RUNX1-sc"]="pink"
  subtypeCol["KMT2A"]="#1F78B5"
  
  subtypeCol["KMT2A_1"]="#1F78B5"
  subtypeCol["KMT2A_2"]="#FF00B6"
  subtypeCol["KMT2A_3"]="#00FF00"
  
  subtypeCol["Ph"]="magenta3"
  subtypeCol["Ph-like"]="red4"
  subtypeCol["Ph-like(CRLF2)"]="red4"
  subtypeCol["Ph-like(non-CRLF2)"]="pink"
  
  subtypeCol["Ph-major"]="magenta3"
  subtypeCol["Ph-minor"]="red4"
  
  subtypeCol["Ph_1"]="magenta3"
  subtypeCol["Ph_2"]="red4"
  
  subtypeCol["PAX5alt"]="#FFA620"
  subtypeCol["PAX5::ETV6"]="#808000"
  subtypeCol["PAX5(P80R)"]="orangered"
  subtypeCol["PAX5 P80R"]="orangered"
  
  subtypeCol["PAX5alt-major"]="#FFA620"
  subtypeCol["PAX5alt-minor"]="pink"
  
  subtypeCol["PAX5alt_1"]="#FFA620"
  subtypeCol["PAX5alt_2"]="pink"
  subtypeCol["PAX5alt_3"]="#5A5156"
  
  subtypeCol["DUX4"]='grey40'
  
  subtypeCol["DUX4_1"]='grey40'
  subtypeCol["DUX4_2"]='#16FF32'
  
  
  subtypeCol["TCF3::PBX1"]="darkgoldenrod4"
  subtypeCol["ZNF384"]="#A8DD00"
  subtypeCol["MEF2D"]="#66C2A6"
  subtypeCol["BCL2/MYC"]="seagreen2"
  subtypeCol["NUTM1"]='black'
  subtypeCol["HLF"]= "skyblue"
  subtypeCol["Hyperdiploid"]="#3E9F32"
  
  subtypeCol["Hyperdiploid_1"]="#3E9F32"
  subtypeCol["Hyperdiploid_2"]="#F6222E"
  
  
  subtypeCol["LowHypo"]="#1E90FF"
  subtypeCol["Low hypodiploid"]="#1E90FF"
  
  subtypeCol["Low hypodiploid_1"]="#1E90FF"
  subtypeCol["Low hypodiploid_2"]="#E4E1E3"
  
  subtypeCol["NearHaploid"]='blue3'
  subtypeCol["Near haploid"]='blue3'
  subtypeCol["Near Haploid"]='blue3'
  
  subtypeCol["iAMP21"]="lightslateblue"
  subtypeCol["IKZF1(N159Y)"]="#CCCC33"
  subtypeCol["IKZF1 N159Y"]="#CCCC33"
  subtypeCol["LowHyper"]="cyan"
  subtypeCol["Bother"]='grey75'
  subtypeCol["Unassigned"]='grey75'
  
  subtypeCol["Low hyperdiploid"]='grey75'
  subtypeCol["Low hyperdiploid_"]='#FF00B6'
  
  subtypeCol["CRLF2(non-Ph-like)"]='grey75'
  subtypeCol["KMT2A-like"]='grey75'
  subtypeCol["ZNF384-like"]='grey75'
  subtypeCol["Other"]='grey75'
  subtypeCol["ZEB2/CEBPE"]="#D27B1C"
  subtypeCol["ZEB2/CEBP"]="#D27B1C"
  
  # subtypeCol["ZEB2/CEBPE"]="orange3"
  # subtypeCol["ZEB2/CEBP"]="orange3"
  
  
  subtypeCol["Y"]="#E6BEFF"
  subtypeCol["CDX2/UBTF"]="#E6BEFF"
  
  # subtypeCol["Unknown"]="#469990"
  subtypeCol["Normal"]="red4"
  subtypeCol["Others"]="grey75"
  # subtypeCol["UNKNOWN"]="grey75"
  # subtypeCol["Unknown"]="grey75"
  
  subtypeCol["UNKNOWN"]="grey75"
  subtypeCol["Unknown"]="grey75"
  
  
  subtypeCol["PAX5alt,Ph-like"]="grey75"
  subtypeCol["Ph-like,iAMP21"]="grey75"
  
  subtypeCol["_Prediction"]="red4"
  
  subtypeCol["AML"]="red4"
  subtypeCol["BALL"]="#FFA620"
  subtypeCol["TALL"]='blue3'
  
  subtypeCol["NormalCD19"]="red4"
  
  subtypeCol["CHLA"]="red4"
  subtypeCol["BALL2239"]='blue3'
  
  subtypeCol["mRNA"]="red4"
  subtypeCol["totalRNA"]='blue3'
  
  
  subtypeCol["TestSample"]="red4"
  
  subtypeCol["No"]="grey90"
  
  subtypeCol["Yes"]="skyblue2"
  subtypeCol["Yes1"]="skyblue2"
  subtypeCol["Yes2"]="tomato2"
  
  # subtypeCol["Female"]="palevioletred1"
  # subtypeCol["Male"]="dodgerblue1"
  #
  # subtypeCol["Adult"]="red4"
  # subtypeCol["Childhood"]="#66C2A6"
  # subtypeCol["NA"]="white"
  
  subtypeCol["European"]="pink"
  subtypeCol["Hispanic"]="#3E9F32"
  subtypeCol["African"]="magenta3"
  subtypeCol["East Asian"]="#1E90FF"
  subtypeCol["South Asian"]="orange"
  subtypeCol["Native Amerian"]="orangered"
  
  subtypeCol["1.goldRef_only"]="#FFA620"
  subtypeCol["2.testing_only"]="#3E9F32"
  subtypeCol["4.notUsed"]="#1E90FF"
  subtypeCol["3.both"]="#E6BEFF"
}

df_ball_col=data.frame(
  diag=names(subtypeCol),
  col=subtypeCol,
  obs=1:length(subtypeCol),
  stringsAsFactors = F
)

Col_PhLikeSubgroup=c()
{
  Col_PhLikeSubgroup["Ph"]="magenta3"
  Col_PhLikeSubgroup["Ph-like"]="red4"
  Col_PhLikeSubgroup["PAX5alt"]="#FFA620"
  Col_PhLikeSubgroup["Hyperdiploid"]="#3E9F32"
  Col_PhLikeSubgroup["Near haploid"]='blue3'
  Col_PhLikeSubgroup["iAMP21"]="lightslateblue"
  Col_PhLikeSubgroup["Low hyperdiploid"]="skyblue"
  Col_PhLikeSubgroup["CRLF2(non-Ph-like)"]="darkgoldenrod4"
  Col_PhLikeSubgroup["Others"]="#1F78B5"
  Col_PhLikeSubgroup["Other"]="#1F78B5"
}

col_list=list(subtypeCol=subtypeCol,Col_PhLikeSubgroup=Col_PhLikeSubgroup)


get_cols=function(in_group,col_default=unique(subtypeCol)){
  cols_in=col_default[1:length(unique(in_group))]
  names(cols_in)=unique(in_group)
  cols_in
}


get_cols_cat=function(value,cols_in,cols_in_default){
  if(!any(names(cols_in) %in% value)){
    print("Supplied color labels not inclued in values levels, use default cols instead")
    if(length(unique(value))>  length(cols_in_default)){message("Defaule color number not enough, the rest will using grey")}
    cols_out=cols_in_default[1:length(unique(value))]
    names(cols_out)=sort(unique(value))
  }
  
  if(any(names(cols_in) %in% value)){
    cols_out=cols_in[names(cols_in) %in% value]
    if(length(cols_out)<length(levels(as.factor(value)))){
      cols_out=cols_in[1:length(levels(as.factor(value)))]
    }
  }
  cols_out
}

get_cols_con=function(value,low,top){}














