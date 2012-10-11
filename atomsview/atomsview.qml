import QtQuick 1.1

 Rectangle {
     id: container
     color: "lightgray"
     
	 property real viewscale: viewHeight / height
     property real viewHeight: 10
     property real viewWidth: width * viewscale
     
	function elementColor(element) {
		var colors = {'H': "white", 'C': "#442211", 'O': "#880000", 'N': "#0000AA"}
		return colors[element]
	}
	
	function elementSize(element) {
		var sizes = {'H': 30, 'C': 45, 'O': 40, 'N':40}
		return sizes[element]
	}
	
	
	Repeater {
		 model: atomsModel
		 
		Rectangle {
			id: atom
			x: (atomx/viewWidth + .5) * container.width - radius
			y: (atomy/viewHeight + .5) * container.height - radius
			z: atomz + 100
			width: radius*2
			height: radius*2
			//radius: covalentRadius/viewscale
			radius: elementSize(element)
			color: elementColor(element)
			border.color: "black"
			border.width: 4
			clip: true
			
			Rectangle {
				id: elementBackground
				visible: type != false
				anchors.centerIn: parent
				color: "#22FFFFFF"
				width: elementName.width + 20
				height: elementName.height + 10
				radius: 15
					
				Text {
					id: elementName
					text: type
					anchors.centerIn: parent
					font.pointSize: 10
				}
			}
			
			MouseArea {
				anchors.fill: parent
				enabled: description || type
				
				onClicked: {
					var dx = mouse.x - (x+radius);
					var dy = mouse.y - (y+radius);
					if (dx^2 + dy^2 < radius^2) {
						tooltipText.text = (description)? description : "No description";
						toolTip.visible = true;
						toolTip.x = atom.x;
						toolTip.y = atom.y;
					}
				}
			}
			
		}
	}
	
	MouseArea {
		id: dragArea
		anchors.fill: parent
		
		property real sensivity: 0.01
		 
		property real lastx: 0
		property real lasty: 0
		 
		onPressed: { 
			lastx = mouse.x; lasty = mouse.y
			toolTip.visible = false
		}
		onPositionChanged: {
			 var dx = mouse.x - lastx;
			 var dy = mouse.y - lasty;
			 lastx = mouse.x;
			 lasty = mouse.y;
			 viewState.rotate(dx*sensivity, dy*sensivity)
		}
	}
	
	Rectangle {
		id: toolTip
		width: tooltipText.width + 20
		height: tooltipText.height + 20
		z: 1000
		radius: 15
		visible: false
		color: "lightgray"
		border.width: 1
		border.color: "black"
		
		Text {
			id: tooltipText
			anchors.centerIn: parent
			text: "tooltip"
		}
		
		MouseArea {
			onPressed: visible = false
		}
	}
 }