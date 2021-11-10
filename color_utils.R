#
# Display several images in a single plot
# by Frédéric Blanchard
#
# Parameters:
#   - image as an array OR filename
#
display<-function(...){
  imageList<-list(...)
  totalWidth<-0
  maxHeight<-0
  for (img in imageList){
    if(is.character(img))
      img<-readPNG(img)
    dimg<-dim(img)
    totalWidth<-totalWidth+dimg[2]
    maxHeight<-max(maxHeight,dimg[1])
  }
  par(mar=c(0,0,0,0))
  plot(c(0,totalWidth),c(0,maxHeight),type="n",asp=1,xaxt="n",yaxt="n",xlab="x",ylab="y")
  offset<-0
  for (img in imageList){
    dimg<-dim(img)
    rasterImage(img,offset,0,offset+dimg[2],dimg[1])
    offset<-offset+dimg[2]
  }
}

#
# Limit pixel values in the [0,1] interval
#
# Parameter:
# - pixel: array containing real numbers
#
# Returns:
#  array containing real numbers between 0 and 1
#
clamp<-function(pixel){
  for(i in 1:length(pixel)){
    if(pixel[i]<0){
      pixel[i]<-0
    }
    else if(pixel[i]>1){
      pixel[i]<-1
    }
  }
  return(pixel)
}

#
# RGB to HSL conversion function
# http://www.rapidtables.com/convert/color/rgb-to-hsl.htm
#
# Parameter:
# - pixel: array containing red, green and blue values
#
# Returns:
#  array containing hue, saturation and lightness values
#
rgb2hsl<-function(pixel){
  r<-pixel[1]
  g<-pixel[2]
  b<-pixel[3]
  
  cmax<-max(r,g,b)
  cmin<-min(r,g,b)
  delta<-cmax-cmin
  
  # Lightness calculation
  l<-(cmax+cmin)/2
  
  # Saturation calculation
  s<-0
  if(delta!=0){
    s<-delta/(1-abs(2*l-1))
  }
  
  # Hue calculation
  h<-0
  if(delta!=0){
    if(cmax==r){
      h<-(((g-b)/delta)%%6)*60
    }
    if(cmax==g){
      h<-(2+(b-r)/delta)*60
    }
    if(cmax==b){
      h<-(4+(r-g)/delta)*60
    }
  }
  return(c(h,s,l))
}

#
# HSL to RGB conversion function
# http://www.rapidtables.com/convert/color/hsl-to-rgb.htm
#
# Parameter:
# - pixel: array containing hue, saturation and lightness values
#
# Returns:
#  array containing red, green and blue values
#
hsl2rgb<-function(pixel){
  h<-pixel[1]
  s<-pixel[2]
  l<-pixel[3]
  
  c<-(1-abs(2*l-1))*s
  x<-c*(1-abs((h/60)%%2-1))
  
  r<-0
  g<-0
  b<-0
  if((h>=0)&(h<60)){
    r<-c
    g<-x
  }
  if((h>=60)&(h<120)){
    r<-x
    g<-c
  }
  if((h>=120)&(h<180)){
    g<-c
    b<-x
  }
  if((h>=180)&(h<240)){
    g<-x
    b<-c
  }
  if((h>=240)&(h<300)){
    r<-x
    b<-c
  }
  if((h>=300)&(h<360)){
    r<-c
    b<-x
  }
  
  m<-l-c/2
  return(clamp(c(r+m,g+m,b+m)))
}
