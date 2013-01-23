import QtQuick 1.1

 Rectangle {
     id: container
     color: "lightgray"
     
	 property real viewscale: viewHeight / height
     property real viewHeight: 10
     property real viewWidth: width * viewscale
     
    function elementColor(element) {
		var col = viewState.jmolColor(element)
		return Qt.rgba(col[0], col[1], col[2], 1)
	}
     
	function elementColor2(element) {
		var colors = {'H': "white", 'C': "#442211", 'O': "#880000", 'N': "#0000BB", 'Al': "darkgray", 'Si': "turquoise"}
		return colors[element]
	}
	
	function elementSize(element) {
		var sizes = {'H': .3, 'C': .45, 'O': .40, 'N': .40, 'Al': .60, 'Mg': .60, 'Si': .45}
		return sizes[element]
	}
	
	
	Repeater {
		model: atomsModel
		
		Item {
			id: atom
			x: (atomx/viewWidth + .5) * container.width
			y: (atomy/viewHeight + .5) * container.height
			z: atomz + 100
		
			Rectangle {
				id: circle
				x: -radius
				y: -radius
				width: radius*2
				height: radius*2
				//radius: covalentRadius/viewscale
				radius: elementSize(element)/viewscale
				color: elementColor(element)
				border.color: "black"
				border.width: 4
				clip: true
				
				Rectangle {
					id: elementBackground
					visible: type != false
					anchors.centerIn: parent
					color: "#33FFFFFF"
					width: elementName.width + 15
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
						var rad = circle.radius;
						var dx = mouse.x - rad;
						var dy = mouse.y - rad;
						if (dx*dx + dy*dy < rad*rad) {
							tooltipText.text = (description)? description : "No description";
							toolTip.visible = true;
							toolTip.x = atom.x;
							toolTip.y = atom.y;
						}
					}
				}
				
			}
		}
	}
	
	
	Repeater {
		model: bondsModel
		
		Item {
			id: bond
			x: ((x1+x2)/2/viewWidth + .5) * container.width
			y: ((y1+y2)/2/viewHeight + .5) * container.height
			z: (z1+z2)/2 + 100
			rotation: Math.atan((y2-y1)/(x2-x1)) *180 / Math.PI
			visible: length3d < 3
			
			property real length2d: hypotenuse(x2-x1, y2-y1)
			property real length3d: hypotenuse(length2d, z2-z1)
			
			property real len1: beamHalfLength(element1)
			property real len2: beamHalfLength(element2)
			
			function hypotenuse(a, b) {
				return  Math.sqrt(a*a + b*b)
			}
			
			function beamHalfLength(element) {
				return length2d * (.5 - elementSize(element)/length3d) / viewscale
			}
			
			Rectangle {
				x: -(x2 > x1? len1 : len2)
				y: -height / 2
				width: len1+len2
				height: 8
				border.color: "black"
				color: "darkgray"
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
			 viewState.rotate(dx*sensivity, -dy*sensivity)
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