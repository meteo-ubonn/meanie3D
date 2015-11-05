if (typeof (M3D) == 'undefined' || M3D == null) {

    M3D = {};

    M3D.graphPercentWidth = 0.85;
    M3D.radius = 12;
    M3D.color = d3.scale.category20();

    /**
     * Class representing a graph of tracks constructed
     * from the data in a tracking dictionary. Also contains
     * helper methods for filtering by id and methods to
     * index or key up the data along criteria.
     */
	function TrackingGraph(dictionary, filterById, filterByTrackLength) {
	    // Data
	    this.dictionary = dictionary;
	    this.selectedId = false;
		this.nodes = dictionary.tree.nodes;
		this.links = dictionary.tree.links;
		// Filtering and indexing
		this.linkedNodes = new Set();
		this.steps = new Array();
		this.ids = new Array();

		this.trackLengthIndex = {
		    values : new Array(),   // distinct track lengths
		    map : new Map()         // ids of the tracks in each class
		}
		this.selectedTrackLength = false;

        /**
         * Sets the id to filter by.
         * @param {number} id
         */
        this.setId = function(filterById) {
            this.selectedId = filterById;
            if (filterById) {
                this.linkedNodes = this.getLinkedNodesByID(parseInt(filterById));
            }
        }

        /**
         * Sets the track length to filter by
         * @param {number} id
         */
        this.setTrackLength = function(length) {
            var self = this;
            this.selectedTrackLength = length;
            if (length) {
                var len = parseInt(this.selectedTrackLength);
                var ids = this.trackLengthIndex.map.get(len);
                if (ids) {
                    // Add all linked nodes of all tracks of the selected length
                    ids.forEach(function(id) {
                        var li = self.getLinkedNodesByID(parseInt(id));
                        li.forEach(function(n) {
                            self.linkedNodes.add(n);
                        });
                    });
                }
            }
        }

        // Filtering

        /** @return a function to filter nodes. */
        this.nodeFilter = function(node) {
            if (this.selectedId || this.selectedTrackLength) {
		        return this.linkedNodes.has(node);
            }
            return true;
        }

        /** @return a function to filter links. */
	    this.linkFilter = function(link) {
	        if (this.selectedId || this.selectedTrackLength) {
                return this.linkedNodes.has(link.source) || this.linkedNodes.has(link.target);
            }
            return true;
        }

        this.filteredNodes = function() {
            var self = this;
            var fn = new Array();
            this.nodes.forEach(function(node) {
                if (self.nodeFilter(node)) {
                    fn.push(node);
                }
            });
            return fn;
        }

        this.filteredLinks = function() {
            var self = this;
            var fl = new Array();
            this.links.forEach(function(link) {
                if (self.linkFilter(link)) {
                    fl.push(link);
                }
            });
            return fl;
        }

        // Indexing and Searching

        /**
         * Turns the uuid references in the links array into
         * object references.
         * @param {array} nodes
         * @param {array} links
         */
        this.resolveLinks = function() {
            var self = this;
            var links = this.links;
            this.nodes.forEach(function(node) {
                self.links.forEach(function(link) {
                    if (link.source == node.uuid) {
                        link.source = node;
                    }
                    if (link.target == node.uuid) {
                        link.target = node;
                    }
                });
            });
        };

        /**
         * Gets all nodes with the matching id
         * @param {Array} nodes
         * @param {number} id
         * @return {Set} of nodes.
         */
        this.nodesWithId = function(id) {
            var nodeSet = new Set();
            this.nodes.forEach(function(node) {
                if (node.id == id)
                    nodeSet.add(node);
            });
            return nodeSet;
        }

        /**
         * Figures out what the id of the node with the given uuid is.
         * @param uuid
         * @returns id
         */
        this.findIdForUuid = function(uuid) {
            var id;
            for (var i=0; i < this.nodes.length && !id; i++) {
                var node = this.nodes[i];
                if (node.uuid == uuid) {
                    id = node.id;
                }
            }
            return id;
        }

        /**
         * Finds a cluster with the given uuid.
         * @param uuid
         * @returns cluster or false.
         */
        this.getClusterForNode = function(node) {
            var cluster = false;
            for (var ti = 0; ti < this.dictionary.tracks.length && !cluster; ti++) {
                var t = this.dictionary.tracks[ti];
                if (t.id == node.id) {
                    for (var ci = 0; ci < t.clusters.length && !cluster; ci++) {
                        var c = t.clusters[ci];
                        if (c.uuid == node.uuid) {
                            cluster = c;
                        }
                    }
                }
            };
            return cluster;
        }

        /**
         * Obtain by id a list of nodes that are linked, starting
         * with all nodes with the given id.
         * @param {Array} nodes
         * @param {Array} links
         * @param {number} id
         * @return {Set} of nodes.
         */
        this.getLinkedNodesByID = function(id) {
            var self = this;
            var roots = this.nodesWithId(id);
            var linked = new Set();
            roots.forEach(function(root) {
                self.getLinkedNodes(root, linked);
            });
            return linked;
        }

        /**
         * Recursively collects all ids of nodes linking
         * to or linked to any node with this id.
         */
        this.getLinkedNodes = function(node,linkedNodes) {
            var self = this;
            if (linkedNodes.has(node)) {
                return;
            }
            linkedNodes.add(node);
            var nextNodes = this.nodesWithId(node.id);
            // Now gather up all nodes that point to or
            // from any of the nodes in the list
            this.links.forEach(function(link) {
                if (nextNodes.has(link.target) || nextNodes.has(link.source)) {
                    var nextNode = null;
                    if (link.target == node) {
                        self.getLinkedNodes(link.source, linkedNodes);
                    } else if (link.source == node) {
                        self.getLinkedNodes(link.target, linkedNodes);
                    }
                }
            });
        }

        /**
         * Indexes the whole tracking dictionary into a few
         * indexes used for navigation.
         */
        this.index = function() {
            var self = this;
            // Index ids/steps
            var ids = this.ids;
            self.nodes.forEach(function(node) {
                if (self.steps.indexOf(node.step) < 0)
                    self.steps.push(node.step);
                if (self.ids.indexOf(node.id) < 0)
                    self.ids.push(node.id);
                });
            this.ids.sort(function(a,b) {
                return parseInt(a.id) < parseInt(b.id);
            });
            this.steps.sort(function(a,b) {
                return parseInt(a.step) < parseInt(b.step);
            });
            // index track dictionary
            this.dictionary.tracks.forEach(function(track) {
                var tli = self.trackLengthIndex.values.indexOf(track.length);
                if ( tli < 0) {
                    self.trackLengthIndex.values.push(track.length);
                }
                var trackIds = self.trackLengthIndex.map.get(track.length);
                if (!trackIds) {
                    trackIds = new Set();
                    self.trackLengthIndex.map.set(track.length,trackIds);
                }
                trackIds.add(track.id);
            });
            // Sort descending
            this.trackLengthIndex.values.sort(function(a,b) {
                return (a > b);
            });
        };

        // Finally index the data
        this.resolveLinks();
		this.index();
		this.setId(filterById);
		if (!filterByTrackLength) {
		    // unless an id was selected, default to the first
		    // entry in the track lengths
		    if (!this.selectedId) {
                this.setTrackLength(this.trackLengthIndex.values[0]);
            }
		} else {
            this.setTrackLength(filterByTrackLength)
        }
    }
    M3D.TrackingGraph = TrackingGraph;


    /**
     * Adds a div#menu to the page, which contains the navigational
     * bubble menu.
     */
    M3D.addBubbleMenu = function(graph) {
        // Extract a linear index as objects
        var children = new Array();
        for (var i = 0; i < graph.trackLengthIndex.values.length; i++) {
            var len = graph.trackLengthIndex.values[i];
            children.push({
                "length": len,
                "count": graph.trackLengthIndex.map.get(len).size
            });
        }
        var dataset = {
            "children": children
        }

        var size = M3D.getMenuSize();
        var diameter = Math.min(size.width, size.height);

        var svg = d3.select("body")
            .select("#menu")
            .append("svg")
            .attr("width", size.width)
            .attr("height", size.height)
            .attr("class", "bubble");

        var format = d3.format(",d");
        var color = d3.scale.category20c();

        var bubble = d3.layout.pack()
            .sort(null)
            .size([diameter, diameter])
            .padding(1.5)
            .value(function(d) { return d.length; })

        var node = svg.selectAll(".node")
            .data(bubble.nodes(dataset))
            .enter().append("g")
            .attr("class", "node")
            .attr("transform", function(d) {
                return "translate(" + d.x + "," + d.y + ")";
            });

        node.append("title")
            .text(function(d) {
                return d.length;
             });

        node.append("circle")
            .attr("r", function(d) { return 10 * Math.log(1+d.length); })
            .style("fill", function(d) { return color(d.length); });

        node.append("text")
            .attr("dy", ".3em")
            .style("text-anchor", "middle")
            .text(function(d) {
                return d.count;
             });

         d3.select(self.frameElement).style("height", diameter + "px");
    }

    /**
     * Adds a div#menu to the page, which contains the navigational
     * bar chart
     */
    M3D.addBarchartMenu = function(graph) {

        // Extract a linear index as objects
        var classes = new Array();
        var maxCount = -1;
        var maxLength = -1;
        for (var i = 0; i < graph.trackLengthIndex.values.length; i++) {
            var len = graph.trackLengthIndex.values[i];
            var count = graph.trackLengthIndex.map.get(len).size;
            classes.push({
                length: len,
                count: count
            });
            if (count > maxCount) maxCount = count;
            if (len > maxLength) maxLength = len;
        }
        classes.sort(function(a,b) {
            return (a.length > b.length);
        });
        var color = d3.scale.ordinal(maxCount);

        var size = M3D.getMenuSize();
        var menu = d3.select("body")
            .select("#menu")
            .select("#barchart")
            .append("div")

        menu.selectAll("div")
            .data(classes)
            .enter()
            .append("div")
            .attr("class", function(d) {
                return d.length == graph.selectedTrackLength ? "bar selected" : "bar";
            })
            .style("width", function(d) {
                var minWidth = 80;
                var width = minWidth + (size.width - minWidth) * (maxCount - d.count) / maxCount - 20;
                return width + "px";
            })
            .style("opacity", function(d) {
                return color(d.count);
            })
            .on("click", function(d) {
                M3D.filterByTrackLength(d.length);
            })
            .text(function(d) {
                return d.length + " (" + d.count + ")";
            })

        if (graph.selectedId) {
            d3.select("#menu")
                .select("#top")
                .append("div")
                .attr("id","trackinfo")
                .text("Selected track #" + graph.selectedId);
        }

    }

    /**
     * Adds a div#graph to the page, which contains the actual graph.
     */
    M3D.addTrackingGraph = function(graph) {
        var size = M3D.getGraphSize();

        var force = d3.layout.force()
            .charge(-200)
            .linkDistance(50)
            .size([size.width, size.height]);

        var svg = d3.select("body")
            .select("#graph")
            .append("svg");

        var nodes = graph.filteredNodes();
        var links = graph.filteredLinks();

        force
            .nodes(nodes)
            .links(links)
            .start();

        // build the arrow.
        svg.append("svg:defs").selectAll("marker")
            .data(["end"])      // Different link/path types can be defined here
            .enter()
            .append("svg:marker")    // This section adds in the arrows
                .attr("id", String)
                .attr("viewBox", "0 -5 10 10")
                .attr("refX", 25)
                .attr("refY", 0)
                .attr("markerWidth", 5)
                .attr("markerHeight", 5)
                .attr("orient", "auto")
                .attr("fill", "#555555")
                .attr("opacity", 0.5)
            .append("svg:path")
                .attr("d", "M0,-5L10,0L0,5");

        var link = svg.append("svg:g").selectAll("path")
            .data(force.links())
            .enter()
            .append("line")
            .attr("class", "link")
            .attr("marker-end", "url(#end)")
            .style("stroke", function(d) {
                if (d.type == 0) {
                    return "lightgrey"
                } else if (d.type == 1) {
                    return "red";
                } else {
                    return "blue";
                }
            })
            .style("stroke-width", function(d) {
                return 2;
            });

        var path = svg.append("svg:g").selectAll("path")
            .data(force.links())
            .enter().append("svg:path")
            .attr("class", "link")
            .attr("marker-end", "url(#end)");

        var node = svg.selectAll(".node")
            .data(nodes)
            .enter().append("circle")
            .attr("class", "node")
            .attr("r", function(d) {
                return 2.5 * Math.log(1.0 + d.size);
            })
            .style("fill", function(d) { return M3D.color(d.id); })
            .on("dblclick", function(d) {
                M3D.filterById(d.id);
            })
            .on("mouseover", function(d) {
                var coordinates = d3.mouse(this);
                M3D.showNodeOverlay(graph,d,coordinates);
            })
            .on("mouseout", function(d) {
                M3D.hideNodeOverlay();
            })
            .call(force.drag);

        var label = svg.selectAll("text")
            .data(nodes)
            .enter().append("text")
            .text(function(d) {
                return d.id;
                // return d.uuid + "/" + d.id;
            })
            .attr("fill","black")
            .attr("text-anchor", "middle")
            .attr("text-align", "center")
            .attr("alignment-baseline","central")
            .attr("font-weight", "bold");

        node.append("title")
            .text(function(d) {
                return d.name;
            });

        force.on("tick", function() {
            link
                .attr("x1", function(d) { return d.source.x; })
                .attr("y1", function(d) { return d.source.y; })
                .attr("x2", function(d) { return d.target.x; })
                .attr("y2", function(d) { return d.target.y; })
                .attr("opacity", function(d) {
                    return 0.25 + 0.75 * (d.source.step / graph.steps.length);
                });


            node
                .attr("cx", function(d) { return d.x; })
                .attr("cy", function(d) { return d.y; })
                .attr("opacity", function(d) {
                    return 0.25 + 0.75 * (d.step / graph.steps.length);
                });

            label
                .attr("x", function(d, i) {
                    return d.x;
                })
                .attr("y", function(d, i) {
                    return d.y;
                });
        });
    }

    M3D.showNodeOverlay = function(graph, node, coordinates) {
        var overlay = d3.select("#node-overlay");
        overlay.attr("display","none");
        var cluster = graph.getClusterForNode(node);
        if (cluster) {
            overlay.select("#node-data").remove();
            var list = overlay
                .append("ul")
                .attr("id","node-data");
            list.append("li").text("id: " + node.id);
            list.append("li").text("uuid: " + cluster.uuid);
            list.append("li").text("size: " + cluster.size);
            list.append("li").text("mins: " + cluster.min);
            list.append("li").text("maxs: " + cluster.max);
            list.append("li").text("medians: " + cluster.median);
            list.append("li").text("margin points: " + cluster.has_margin_points);
            list.append("li").text("source: " + cluster.sourcefile);
            list.append("li").text("center: " + cluster.geometrical_center);
        }

        var x = coordinates[0] + 200;
        var y = coordinates[1] - 200;

        overlay
            .style("background-color", M3D.color(node.id))
            .style("display","block")
            .style("opacity", 0.8)
            .style("position","absolute")
            .style("left",x + "px")
            .style("top", y + "px");
    }

    M3D.hideNodeOverlay = function() {
        var overlay = d3.select("#node-overlay");
        overlay.style("display","none");
    }


    M3D.getViewportSize = function() {
        var size = {
            'width' : Math.max(document.documentElement.clientWidth, window.innerWidth || 0),
            'height' : Math.max(document.documentElement.clientHeight, window.innerHeight || 0)
        }
        return size;
    }

    M3D.getGraphSize = function() {
        var viewport = M3D.getViewportSize();
        var size = {
            'width' : M3D.graphPercentWidth * viewport.width,
            'height' : viewport.height
        }
        return size;
    }

    M3D.getMenuSize = function() {
        var viewport = M3D.getViewportSize();
        var size = {
            'width' : (1.0 - M3D.graphPercentWidth) * viewport.width,
            'height' : viewport.height
        }
        return size;
    }

    /**
     * Courtesy of https://css-tricks.com/snippets/javascript/get-url-variables/
     */
    M3D.getUrlParameter= function(param)
    {
       var query = window.location.search.substring(1);
       var vars = query.split("&");
       for (var i=0;i<vars.length;i++) {
               var pair = vars[i].split("=");
               if(pair[0] == param){return pair[1];}
       }
       return(false);
    }

    M3D.resize = function() {
        var menuSize = M3D.getMenuSize();
        d3.select("#menu")
            .style("width", menuSize.width)
            .style("height", menuSize.height);

        var graphSize = M3D.getGraphSize();

        d3.select("#graph").select("svg")
            .attr("width", graphSize.width)
            .attr("height", graphSize.height);

    }

    /**
     * Reloads the current page with an id for filtering.
     * @param A node or link or any other object with property 'id'.
     */
    M3D.filterById = function(d) {
        var location = window.location.toString();
        var i = location.indexOf("?");
        if ( i > 0) {
            location = location.substring(0,i);
        }
        window.location.href = location + "?id="+d;
    }

    M3D.filterByTrackLength = function(d) {
        var location = window.location.toString();
        var i = location.indexOf("?");
        if ( i > 0) {
            location = location.substring(0,i);
        }
        window.location.href = (location + "?tracklength=" + d);
    }

    M3D.search = function(graph,id,uuid) {
        if (id) {
            M3D.filterById(id);
        } else if (uuid) {
            var _id = graph.findIdForUuid(uuid);
            M3D.filterById(_id);
        }
        return false;
    }

}