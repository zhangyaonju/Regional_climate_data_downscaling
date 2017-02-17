tile_climate_data_ncep<-function(year,h,v){
  
  parameter.file<-paste('/data/ifs/users/yzhang/NCEP/h',
                        formatC(h,format='d',width=2,flag='0'),'v',formatC(v,format='d',width=2,flag='0'),'.tif',sep='')
  
  input.8day.file<-c(paste('/data/ifs/users/yzhang/NCEP/temperature/air.2m.8day.',year,'.nc',sep=''),
                     paste('/data/ifs/users/yzhang/NCEP/radiation/rad.8day.',year,'.nc',sep=''))
  
  air.temp.raw.data<-nc_open(input.8day.file[1])
  air.temp.8day.data<-ncvar_get(air.temp.raw.data,'air')
  
  
  radiation.raw.data<-nc_open(input.8day.file[2])
  radiation.8day.data<-ncvar_get(radiation.raw.data,'dswrf')
   
  
  resample.parameter.data<-readGDAL(parameter.file)
  coor.orig<-summary(resample.parameter.data)[[5]]
  resample.parameter<-resample.parameter.data@data
  coor.x<-array(NA,dim=c(2400))
  coor.y<-array(NA,dim=c(2400))
  for(i in 1:2400){
    coor.x[i]<-coor.orig[1,1]+(i-1)*463.3127
    coor.y[i]<-coor.orig[2,1]+(2400-i)*463.3127
  }
  #proj<-array(,dim=c(2400*2400,2))
  #for(i in 1:2400){
  #  for(j in 1:2400){
  #    proj[2400*(i-1)+j,1]<-coor.x[i]
  #    proj[2400*(i-1)+j,2]<-coor.y[j]
  #  }
  #}
  #reproject<-project(proj,summary(resample.parameter.data)$proj4string,inv=TRUE)
  
  
  radiation.fine.output<-paste('/data/eomf/users/yzhang/NCEP/tile/radiation/',year,'/rad.h',
                               formatC(h,format='d',width=2,flag='0'),'v',formatC(v,format='d',width=2,flag='0'),'.',
                               year,'.nc',sep='')
  air.temp.fine.output<-paste('/data/eomf/users/yzhang/NCEP/tile/temperature/',year,'/temp.h',
                              formatC(h,format='d',width=2,flag='0'),'v',formatC(v,format='d',width=2,flag='0'),'.',
                              year,'.nc',sep='')
  
  fine.air.temp.data<-array(,dim=c(2400,2400,46))
  fine.radiation.data<-array(,dim=c(2400,2400,46))
  
  for (scene.num in 1:46){
    fine.air.temp.data[,,scene.num]<-air.temp.8day.data[,,scene.num][resample.parameter[,5]+193]*resample.parameter[,1]+air.temp.8day.data[,,scene.num][resample.parameter[,5]+194]*resample.parameter[,2]+
      air.temp.8day.data[,,scene.num][resample.parameter[,5]+386]*resample.parameter[,3]+air.temp.8day.data[,,scene.num][resample.parameter[,5]+385]*resample.parameter[,4]
    fine.radiation.data[,,scene.num]<-radiation.8day.data[,,scene.num][resample.parameter[,5]+193]*resample.parameter[,1]+radiation.8day.data[,,scene.num][resample.parameter[,5]+194]*resample.parameter[,2]+
      radiation.8day.data[,,scene.num][resample.parameter[,5]+386]*resample.parameter[,3]+radiation.8day.data[,,scene.num][resample.parameter[,5]+385]*resample.parameter[,4]
  }
  
  
  x<-ncdim_def('eastward distance from southwest corner of domain in projection coordinates','m',coor.x)
  y<-ncdim_def('northward distance from southwest corner of domain in projection coordinates','m',coor.y)

  t<-ncdim_def("Time",'hours since 1800-01-01 00:00:0.0',radiation.raw.data$dim[[3]]$vals,unlim=TRUE)
  

  ncrad<-ncvar_def('dswrf',"W/m^2",list(x,y,t),missval=-9.96920996838687e+36,longname='8-day Downward Solar Radiation Flux at surface',prec='float')

  
  ncrad.out<-nc_create(radiation.fine.output,ncrad)

  ncvar_put(ncrad.out,'dswrf',fine.radiation.data)

  nc_close(ncrad.out)
  
  
  ####=================================================================
  

  
  x<-ncdim_def('eastward distance from southwest corner of domain in projection coordinates','m',coor.x)
  y<-ncdim_def('noethward distance from southwest corner of domain in projection coordinates','m',coor.y)
  t<-ncdim_def("Time",'hours since 1800-01-01 00:00:0.0',air.temp.raw.data$dim[[3]]$vals,unlim=TRUE)
  
  ncair<-ncvar_def('air',"degK",list(x,y,t),missval=-9.96920996838687e+36,longname='8-day 2m air temperature',prec='float')

  
  ncair.out<-nc_create(air.temp.fine.output,ncair)
  
  ncvar_put(ncair.out,'air',fine.air.temp.data)

  
  nc_close(ncair.out)

}