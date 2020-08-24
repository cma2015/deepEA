if(window.attachEvent){
		window.attachEvent("onload", test)
	} else {
		window.addEventListener("load", test)
	}
		function test(){
			function $(param){
				if(arguments[1] == true){
					return document.querySelectorAll(param);
				}else{
					return document.querySelector(param);
				}
			}
			var $box = $(".box");
			var $box1 = $(".box-1 ul li",true);
			var $box2 = $(".box-2 ul");
			var $box3 = $(".box-3");
			var $length = $box1.length;
			
			var str = "";
			for(var i =0;i<$length;i++){
				if(i==0){
					str +="<li class='on'>"+(i+1)+"</li>";
				}else{
					str += "<li>"+(i+1)+"</li>";
				}
			}
			$box2.innerHTML = str;
			
			var current = 0;
			
			var timer;
			timer = setInterval(go,3000);
			function go(){
				for(var j =0;j<$length;j++){
					$box1[j].style.display = "none";
					$box2.children[j].className = "";
				}
				if($length == current){
					current = 0;
				}
				$box1[current].style.display = "block";
				$box2.children[current].className = "on";
				current++;
			}
			
			for(var k=0;k<$length;k++){
				$box1[k].onmouseover = function(){
					clearInterval(timer);
				}
				$box1[k].onmouseout = function(){
					timer = setInterval(go,1000);
				}
			}
			for(var p=0;p<$box3.children.length;p++){
				$box3.children[p].onmouseover = function(){
					clearInterval(timer);
				};
				$box3.children[p].onmouseout = function(){
					timer = setInterval(go,1000);
				}
			}
			
			for(var u =0;u<$length;u++){
				$box2.children[u].index  = u;
				$box2.children[u].onmouseover = function(){
					clearInterval(timer);
					for(var j=0;j<$length;j++){
						$box1[j].style.display = "none";
						$box2.children[j].className = "";
					}
					this.className = "on";
					$box1[this.index].style.display = "block";
					current = this.index +1;
				}
				$box2.children[u].onmouseout = function(){
					timer = setInterval(go,1000);
				}
			}
			
			$box3.children[0].onclick = function(){
				back();
			}
			$box3.children[1].onclick = function(){
				go();
			}
			function back(){
				for(var j =0;j<$length;j++){
					$box1[j].style.display = "none";
					$box2.children[j].className = "";
				}
				if(current == 0){
					current = $length;
				}
				$box1[current-1].style.display = "block";
				$box2.children[current-1].className = "on";
				current--;
			}
		}