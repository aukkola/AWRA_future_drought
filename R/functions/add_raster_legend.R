add_raster_legend2 <-function(cols,limits,labelss=NULL,x='bottomleft',add=TRUE,
                              spt.cex=2,pch=15,main_title='', title.cex=1.8,
                              plot_loc=c(0.22,0.78,0.07,0.10), clip=FALSE, ysp_title_old=FALSE, ...) {
  
  cal_frac_of_plot_covered <- function(index=1:2) {
    xlims=par("usr")[index]
    xlims[3]=xlims[2]-xlims[1]
    return(xlims)
  }
  
  if (add) {
    xlims=cal_frac_of_plot_covered(1:2)
    ylims=cal_frac_of_plot_covered(3:4)
  } else {
    xlims=c(0,1,1)
    ylims=c(0,1,1)
    plot(c(0,1),c(0,1))
  }
  nlims=length(limits)+1
  
  x=xlims[1]+plot_loc[1]*xlims[3]+(plot_loc[2]-plot_loc[1])*xlims[3]*(matrix(1:nlims,nlims,1)-0.5)/nlims
  
  y=c(ylims[1]+ylims[3]*plot_loc[3],ylims[1]+ylims[3]*plot_loc[4])
  
  marTemp=par("mar")
  par(mar=c(0,0,0,0), xpd=TRUE)
  
  image(x=x,y=y,z=matrix(1:nlims,nlims,2),col=cols[1:nlims], add=TRUE,
        xaxt='n',yaxt='n',xlab='',ylab='', useRaster=TRUE) ###xpd=TRUE
  
  #If legend needs to be outside of plotting region, xpd=NA does not
  #allow this for some reason when using image. This is a hacky solution
  #https://stackoverflow.com/questions/50558656/how-to-plot-image-outside-plot-area-with-parxpd-t
  
  if (clip) {
    
    grid::grid.clip()
    image(x=x,y=y,z=matrix(1:nlims,nlims,2),col=cols[1:nlims], add=TRUE,
          xaxt='n',yaxt='n',xlab='',ylab='', useRaster=TRUE) ###xpd=TRUE
  }
  
  par(mar=marTemp)  
  
  
  dx=(x[2]-x[1])/2
  xp=c(min(x)-dx,min(x)-dx,max(x)+dx,max(x)+dx)
  
  yp=c(y[1]-diff(y)/2,y[2]+diff(y)/2)
  polygon(xp,c(yp,rev(yp)),xpd=NA) ##xpd=TRUE
  #===============================
  ytop=yp[2]+diff(y)
  ybottom=yp[1]-diff(y)
  xt=c(x[1]-dx,x[1],x+dx)
  yt=c(ytop,rep(c(ytop,ybottom),length.out=length(xt)-1))
  if (is.null(labelss)) {
    labelss=c("","",limits,"")
    #  if (limits[1]==0) {
    #   labelss[1]=""
    #  labelss[2]='0'
    # labelss[3]=""
    #}
  } else {
    if(class(labelss)=='list' && length(labelss)==1) labelss=labelss[[1]]
    
    if (limits[1]==0) {
      labelss=c("",labelss[1],"",labelss[-1],"")
    } else {
      labelss=c(labelss[1],"",labelss[-1],"")
    }
    
    if (length(labelss)==(length(xt)-2)) labelss=c("",labelss,"")
    if (length(labelss)==(length(xt)-1)) labelss=c("","",labelss[-2])
  }
  if (length(labelss)==0) labelss=rep("",length(xt))
  text(x=xt,y=yt,labelss,xpd=NA, cex=spt.cex)
  #changed y=max(yt+max(yt*0.15)) from y=max(yt-max(yt*0.35))
  #Had to change this for some plots, but change doesn't work for some others...
 # if (ysp_title_old) {
    text(x=mean(xt[2:length(xt)]),y=1.9, main_title, xpd=NA, cex=title.cex)
#  } else {
 #   text(x=mean(xt),y=max(yt+max(yt*0.1)),main_title,xpd=NA, cex=title.cex)
#  }
  #===============================
}


