
filter.html=function(html){
	html=html[!html==""]
	name=grep("href=/cgi-bin/", html)
	if(length(name)){
		head=max(name)+1
		tail=grep("Citation: ", html)
		if(!length(tail)) tail=length(html)
		return(html[head:tail])
	} else {
		return(NULL)
	}
}

htmlstring=function(html){
	unlist(strsplit(as(html,"character"),"\n"))
}
