import sys
import math
import glob
import ezdxf
import numpy
import logging
import matplotlib
from matplotlib.path import Path
from matplotlib.transforms import Affine2D
from ezdxf.math import Vec2, ConstructionArc
from functools import partial

# used file specs:
# gerber - http://www.ucamco.com/Portals/0/Documents/Ucamco/RS-274X_Extended_Gerber_Format_Specification_201201.pdf
# dxf 2012 - http://images.autodesk.com/adsk/files/acad_dxf1.pdf

default_thickness = 0.01
autofill = False
epsilon = sys.float_info.epsilon


def perpendicular( a ) :
    return Vec2((a[1],-a[0]))


def intersection(s1,e1,s2,e2):
    d = (s1[0]-e1[0])*(s2[1]-e2[1]) - (s1[1]-e1[1])*(s2[0]-e2[0])
    if math.fabs(d) < 2*epsilon: return None
    
    f1 = s1[0]*e1[1]-s1[1]*e1[0]
    f2 = s2[0]*e2[1]-s2[1]*e2[0]
    
    return ((s2-e2)*f1-(s1-e1)*f2)/d


def convert_boundery_to_path(points):
    vertices = list()
    codes = list()
    for segment in points:
        ptype,start,end = segment[:3]

        if ptype == "straight":#ray line intersection with infinite ray from (p_x,p_y) to (inf,p_y) of first point of A
            v,c = next(Path((start,end)).iter_segments())
            vertices.append((v[0],v[1]))
            codes.append(c)
        elif ptype == "arc":#ray arc intersection with infinite ray from (p_x,p_y) to (p_x+1,p_y) of first point of A
            rel,ccw,arc = segment[3:]
            trans   = Affine2D().scale(rel.magnitude).translate(arc.center.x, arc.center.y)
            path    = Path.arc(arc.start_angle, arc.end_angle)# if ccw else Path.arc(arc.end_angle, arc.start_angle,)
            path    = path.transformed(trans)

            for v,c in path.iter_segments(curves=False):
                vertices.append((v[0],v[1]))
                codes.append(c)
        else:
            logging.warning("unknown boundary type found (convert_boundery_to_path)")
    codes = [1] + [2]*(len(vertices)-1)
    ret = Path(vertices,codes=codes, closed=True)#.cleaned(simplify=True)#     
    return ret


#######################################################################
# converting + translation
#######################################################################

def convert_circle(e):
    layer   = e.dxf.layer
    x,y     = e.dxf.center[:2]
    r       = e.dxf.radius*2
    device  = "C,%f*"%r if autofill else "C,%fX%f*"%(r,r-default_thickness)
    
    return [{"layer":layer,"device":device,"x":x,"y":y,"style":"D03"}],True


def convert_arc(e):
    layer               = e.dxf.layer
    center_x,center_y   = e.dxf.center[:2]
    thickness           = e.dxf.thickness
    radius              = e.dxf.radius
    start_angle         = e.dxf.start_angle/180.0*math.pi
    end_angle           = e.dxf.end_angle/180.0*math.pi
    
    rel_x       = radius * math.cos(start_angle)
    rel_y       = radius * math.sin(start_angle)
    start_x     = center_x + rel_x
    start_y     = center_y + rel_y
    end_x       = center_x + radius * math.cos(end_angle)
    end_y       = center_y + radius * math.sin(end_angle)
    
    objects = []
    device = "C,%f*"%(thickness if thickness else default_thickness)
    objects.append({"layer":layer,"command":"G75*"})
    objects.append({"layer":layer,"device":device,"mode":"G01","x":start_x,"y":start_y,"style":"D02"})
    objects.append({"layer":layer,"device":device,"mode":"G03","style":"D01","x":end_x,"y":end_y,"i":-rel_x,"j":-rel_y})
    objects.append({"layer":layer,"command":"G01*"})

    return objects,True


def convert_ellipse(e):
    layer               = e.dxf.layer
    thickness           = e.dxf.thickness
    center_x,center_y   = e.dxf.center[:2]
    major_x,major_y     = e.dxf.major_axis[:2]
    minor_to_major      = e.dxf.ratio
    start_angle         = e.dxf.start_param
    end_angle           = e.dxf.end_param

    tilt        = math.atan(major_y/major_x) if major_x else math.copysign(math.pi/2,major_y)
    sampling    = int((end_angle - start_angle)*8/math.pi)#a full cycle will be sampled with 2pi/(pi/8)= 16 polylines
    device      = "C,%f*"%(thickness if thickness else default_thickness)
    major_len   = math.sqrt(major_x**2+major_y**2)
    minor_len   = (major_x**2+major_y**2) * (minor_to_major**2)
    objects     = []
    
    if autofill:
        objects.append({"layer":layer,"command":"G36*"})
    
    for step in range(sampling):
        t = (end_angle - start_angle) * step / sampling
        x = center_x + major_len*math.cos(t)*math.cos(tilt) - minor_len*math.sin(t)*math.sin(tilt)
        y = center_y + major_len*math.cos(t)*math.sin(tilt) - minor_len*math.sin(t)*math.cos(tilt)
        objects.append({"layer":layer,"device":device,"x":x,"y":y,"style":"D01"})
        
    objects[0]["style"] = "D02"
    
    if autofill:
        objects.append({"layer":layer,"command":"G37*"})
       
    return objects,True


def convert_line(e):
    layer           = e.dxf.layer
    thickness       = e.dxf.thickness
    start_x,start_y = e.dxf.start[:2]
    stop_x,stop_y   = e.dxf.end[:2]
    
    device = "C,%f*"%(thickness if thickness else default_thickness)
    return [{"layer":layer,"device":device,"x":start_x,"y":start_y,"style":"D02"},
            {"layer":layer,"device":device,"x":stop_x,"y":stop_y,"style":"D01"}],True


def convert_insert(e, converters, converted_blocks=dict()):
    block_name      = e.dxf.name
    layer           = e.dxf.layer
    xoffset,yoffset = e.dxf.insert[:2]
    xscale          = e.dxf.xscale
    yscale          = e.dxf.yscale
    rotation        = e.dxf.rotation/180*math.pi
    switch_ccw      = xscale*yscale < 0

    objects = list()
    
    if block_name not in converted_blocks:
        converted_blocks[block_name] = convert_block(e.doc.blocks[block_name], converters)

    for obj in converted_blocks[block_name]:
        obj = dict(obj)#create a copy so we can put the block on multiple places
        if "x" in obj:
            x = xscale*obj["x"]
            y = yscale*obj["y"]
            obj["x"] = (x*math.cos(rotation) - y*math.sin(rotation)) + xoffset
            obj["y"] = (x*math.sin(rotation) + y*math.cos(rotation)) + yoffset
        if "i" in obj:
            i = xscale*obj["i"]
            j = yscale*obj["j"]
            obj["i"] = (i*math.cos(rotation) - j*math.sin(rotation))
            obj["j"] = (i*math.sin(rotation) + j*math.cos(rotation))
        
        if block_name[0] != "*":
            obj["layer"] = layer
        
        if switch_ccw and "mode" in obj:
            if   obj["mode"] == "G03": obj["mode"] = "G02"
            elif obj["mode"] == "G02": obj["mode"] = "G03"
        objects.append(obj)                
    return objects,True


def convert_lwpolyline(e):
    objects = []
    closed  = e.closed
    layer   = e.dxf.layer
    
    for index in range(e.dxf.count - (not e.closed)):
        sx, sy, start_width, end_width, start_bulge = e[index]
        ex, ey, _, _, end_bulge                     = e[(index+1)%e.dxf.count]

        if e.dxf.const_width:
            start_width = e.dxf.const_width
            end_width   = e.dxf.const_width

        start               = Vec2((sx,sy))
        end                 = Vec2((ex,ey))
        contours            = start_width!=0 or end_width!=0
        device              = ("C,%f*"%start_width) if start_width else "C,%f*"%default_thickness
        ortholine           = perpendicular(end - start)
        ortholine_l         = ortholine.magnitude
        
        if start_bulge:
            #rel = ezdxf.math.bulge_center(start, end, bulge)-ezdxf.math.Vec2(start)
            opening_angle   = math.atan(float(start_bulge))*4
        
            rel_start   = (end - start)/2 - ortholine / (2 * math.tan(opening_angle/2))
            rel_end     = start - end + rel_start

            radius      = rel_start.magnitude
            
            start1      = start - rel_start * 0.5 * start_width / radius
            start2      = start + rel_start * 0.5 * start_width / radius
            end1        = end - rel_end * 0.5 * end_width / radius
            end2        = end + rel_end * 0.5 * end_width / radius
        elif ortholine_l:
            start1      = start - ortholine * 0.5 * start_width / ortholine_l
            start2      = start + ortholine * 0.5 * start_width / ortholine_l
            end1        = end - ortholine * 0.5 * end_width / ortholine_l
            end2        = end + ortholine * 0.5 * end_width / ortholine_l

            #to ensure that no gaps are created between segments we need to
            #find the intersection points with the previous segment
            #NOTE: currently disable that part
##            if w1 and w2 and i > 0 and not points[i-1][1]:
##                start3,end3,start4,end4 = lines[-1][:4]
##                ret1 = intersection(start1,end1,start3,end3)
##                ret2 = intersection(start2,end2,start4,end4)
##                
##                if ret1 != None: 
##                    start1 = ret1
##                    lines[-1][1] = ret1
##                if ret2 != None: 
##                    start2 = ret2
##                    lines[-1][3] = ret2
        else:
            start1      = start
            start2      = start
            end1        = end
            end2        = end
        
        if start_bulge:
            opening_angle   = math.atan(float(start_bulge))*4
            mode            = "G03" if float(start_bulge) > 0 else "G02"
            inv_mode        = "G02" if float(start_bulge) > 0 else "G03"
            
            objects.append({"layer":layer,"command":"G75*"})
            if contours:
                ortho1      = perpendicular(end1 - start1)
                ortho2      = perpendicular(end2 - start2)
                cycle1      = -(end1 - start1)/2 - ortho1 / (2 * math.tan(opening_angle/2))
                cycle2      = (end2 - start2)/2 - ortho2 / (2 * math.tan(opening_angle/2))
                
                objects.append({"layer":layer,"command":"G36*"})
                objects.append({"layer":layer,"device":"C,0.001*","mode":"G01","x":start1[0],"y":start1[1],"style":"D02"})
                objects.append({"layer":layer,"device":"C,0.001*","mode":"G01","x":start2[0],"y":start2[1],"style":"D01"})
                objects.append({"layer":layer,"device":"C,0.001*","mode":mode,"style":"D01","x":end2[0],"y":end2[1],"i":cycle2[0],"j":cycle2[1]})
                objects.append({"layer":layer,"device":"C,0.001*","mode":"G01","x":end1[0],  "y":end1[1],  "style":"D01"})
                objects.append({"layer":layer,"device":"C,0.001*","mode":inv_mode,"style":"D01","x":start1[0],"y":start1[1],"i":cycle1[0],"j":cycle1[1]})
                objects.append({"layer":layer,"command":"G37*"})
            else:#circular arc with little round ends but who cares?
                rel = -(start - end)/2 - ortholine / (2 * math.tan(opening_angle/2))
                objects.append({"layer":layer,"device":device,"mode":"G01","x":start[0],"y":start[1],"style":"D02"})
                objects.append({"layer":layer,"device":device,"mode":mode,"style":"D01","x":end[0],"y":end[1],"i":rel[0],"j":rel[1]})
            objects.append({"layer":layer,"command":"G01*"})
        elif ortholine_l > 0:
            if contours:#using a rectangluar approach
                objects.append({"layer":layer,"command":"G36*"})
                objects.append({"layer":layer,"device":"C,0.001*","x":start1[0],"y":start1[1],"style":"D02"})
                objects.append({"layer":layer,"device":"C,0.001*","x":start2[0],"y":start2[1],"style":"D01"})
                objects.append({"layer":layer,"device":"C,0.001*","x":end2[0],  "y":end2[1],  "style":"D01"})
                objects.append({"layer":layer,"device":"C,0.001*","x":end1[0],  "y":end1[1],  "style":"D01"})
                objects.append({"layer":layer,"device":"C,0.001*","x":start1[0],"y":start1[1],"style":"D01"})
                objects.append({"layer":layer,"command":"G37*"})
            else:#straight line with little round ends but who cares?
                if not ortholine_l:
                    print("warning: found a thick polygon line with a length of %d. cannot handle that correctly=>fallback to simple line."%ortholine_l)
                    #print i,count,points[i],points[i+1]
                objects.append({"layer":layer,"device":device,"x":start[0],"y":start[1],"style":"D02"})
                objects.append({"layer":layer,"device":device,"x":end[0],  "y":end[1],  "style":"D01"})
                
    #NOTE: try really closed poly and remove doubled lines
    if autofill and closed:
        objects = [o for o in objects if o.get("command") not in ("G36*","G37*")]
        first = True
        for o in objects:
            if "style" in o:
                o["style"] = "D02" if first else "D01"
                first = False
        objects.insert(0,{"layer":layer,"command":"G36*"})
        objects.append({"layer":layer,"command":"G37*"})

    return objects, True


def convert_boundary(layer, points, fill=False):
    last    = None
    ret     = list()
    if fill: ret.append({"layer":layer,"command":"G36*"})
    for p in points:
        typ,start,end = p[:3]
        if p[0] == "straight":
            if not last or not ezdxf.math.is_close_points(last, start):
                ret.append({"layer":layer,"device":"C,0.001*","x":start.x,"y":start.y,"style":"D02"})
            ret.append({"layer":layer,"device":"C,0.001*","x":end.x,"y":end.y,"style":"D01"})
        else:#arc
            rel,ccw,arc = p[3:]
            
            ret.append({"layer":layer,"command":"G75*"})
            if not last or not ezdxf.math.is_close_points(last, start):
                ret.append({"layer":layer,"device":"C,0.001*","mode":"G01","x":start.x,"y":start.y,"style":"D02"})
            ret.append({"layer":layer,"device":"C,0.001*","mode":"G03" if ccw else "G02","style":"D01","x":end.x,"y":end.y,"i":rel.x,"j":rel.y})
            ret.append({"layer":layer,"command":"G01*"})
        last = end
    if fill: ret.append({"layer":layer,"command":"G37*"})
    return ret

    
def convert_hatch(e):
    pattern_name    = e.dxf.pattern_name
    layer           = e.dxf.layer
    solid           = e.dxf.solid_fill
    hatchstyle      = e.dxf.hatch_style

    boundaries      = list()
    roots           = list()
    ret             = list()

    logging.debug(f"hatchstyle {hatchstyle}")#hatchstyle is currently ignored

    #darkfield = True
    for no,part in enumerate(e.paths.paths):
        external    = part.path_type_flags & 1
        outermost   = part.path_type_flags & 16
        points      = list()
        
        if type(part) is ezdxf.entities.PolylinePath:
            is_closed = part.is_closed
            has_bulge = len(part.vertices[0]) == 3
            
            for index in range(len(part.vertices) - (not is_closed)):
                bulge   = part.vertices[index][2]
                start   = Vec2(part.vertices[index][:2])
                end     = Vec2(part.vertices[(index+1)%len(part.vertices)][:2])
                
                logging.debug(f"polyline start-end {start} {end} {end} {bulge}")
                if bulge:
                    radius = ezdxf.math.bulge_radius(start, end, bulge)
                    arc = ezdxf.math.ConstructionArc.from_2p_radius(start, end, radius,bulge > 0)
                    rel = arc.center - start
                    #print(f"bulge ({start.x:6.2f}, {start.y:6.2f}) ({end.x:6.2f}, {end.y:6.2f}) ({rel.x:6.2f}, {rel.y:6.2f})",bulge)
                    points.append(("arc",start,end,rel,bulge > 0, arc))
                else:
                    points.append(("straight",start, end))
        else:
            for e in part.edges:
                if e.EDGE_TYPE == "LineEdge":
                    points.append(("straight",Vec2(e.start),Vec2(e.end)))
                    
                elif e.EDGE_TYPE == "ArcEdge":
                    arc = ezdxf.math.ConstructionArc(e.center, e.radius, e.start_angle, e.end_angle, e.is_counter_clockwise)
                    rel = arc.center - arc.start_point
                    print("rel")
                    points.append(("arc",arc.start_point,arc.end_point,rel,e.is_counter_clockwise, arc))
                else:
                    logging.warning(f"unhandled edgetype {edgetype}. Try to ignore it but mostelike it will fail.")
            
        if not ezdxf.math.is_close_points(points[0][1],points[-1][2]):
            logging.warning(f"not closed objects are NOT supported within a hatch. distance check failed. Adding straight line to close the gap!")
            points.append(("straight",points[-1][2],points[0][1]))

        #assert(not external or not outermost)
        path = convert_boundery_to_path(points)
        if path is None: continue
        if external:    roots.append((points, path, list()))
        else:           boundaries.append((points, path))
        
    #building dependency/lies within tree
    if not len(roots):
        logging.warning("empty hatch found.")
        return [],True

    while len(boundaries):
        branch      = roots
        points,path = boundaries.pop()
        searching   = True
        children    = list()
        
        while searching:
            searching = False
            for leaf in branch[:]:
                otherpoints,otherpath,subobj = leaf
                if path.contains_point(otherpoints[0][1]):#other lies within our boundary
                    branch.remove(leaf)
                    children.append(leaf)
                elif otherpath.contains_point(points[0][1]):#our boundary lies within other
                    branch      = subobj
                    searching   = True
                    break
            if not searching:
                branch.append((points, path, children))

                

    darkfield   = True
    drawlist    = roots
    while len(drawlist):
        new_drawlist = list()

        for points, path, subobjs in drawlist:
            new_drawlist.extend(subobjs)
        
            last        = None
            last_parity = darkfield

            #debug
##            segments = path.to_polygons()
##            points = list()
##            for s in segments:
##                start = s[-1]
##                for p in s:
##                    points.append(("straight", Vec2(start), Vec2(p)))
##                    start = p
            
            ret.extend(convert_boundary(layer,points,True))
            
        darkfield   = not darkfield
        drawlist    = new_drawlist
        ret.append({"layer":layer,"command":"%LPD*%" if darkfield else "%LPC*%"})
        if not darkfield:
            logging.warning("Using LPC (clear area) command which can potential erase other structures as well. Future version has to replace that by unifing the bounderies to one with additional lines")
        
    ret.append({"layer":layer,"command":"%LPD*%"})
    return ret,False


def convert_block(block, converter):  
    logging.debug(f"converting block {block.name}")
    block_objects = []
    for element in block:
        if element.dxftype() == "HATCH":
            logging.warning("found a HATCH within a INSERT which can cause malfunctional drawing (elements can be erased)")

        child_objects,append = converter[element.dxftype()](element)
            
        if append:  block_objects = block_objects + child_objects
        else:       block_objects = child_objects + block_objects

    return block_objects


def write_gerber(layers):
    for l,objects in layers.items():
        devices = {obj["device"]:None for obj in objects if "device" in obj}
        with open(l+".gbr","w") as file:
            file.write("%FSLAX36Y36*%\n")
            file.write("%MOMM*%\n")
            file.write("%LN"+ l + "*%\n")
        
            counter = 10
            for d in devices.keys():
                file.write("%ADD" + str(counter) + d + "%\n")
                devices[d] = counter
                counter += 1
                
            current_device = None
            for obj in objects:
                if "command" in obj:
                    file.write(obj["command"] +"\n")
                else:
                    if obj["device"] != current_device:
                        file.write("G54D%d*\n"%devices[obj["device"]])
                        current_device = obj["device"]
                        
                    line = ""
                    
                    if "mode" in obj:   line += obj["mode"]
                    if "x" in obj:      line += "X"+str(int(obj["x"] * 1000000))
                    if "y" in obj:      line += "Y"+str(int(obj["y"] * 1000000))
                    if "i" in obj:      line += "I"+str(int(obj["i"] * 1000000))
                    if "j" in obj:      line += "J"+str(int(obj["j"] * 1000000))
                    if "style" in obj:  line += obj["style"]
                    
                    file.write(line+"*\n")
            file.write("M02*\n")
    print("done!")
            

def convert_dxf2gerber(filename):
    dxf = ezdxf.readfile(filename)
    
    converters = {  'CIRCLE'    : convert_circle,
                    'ARC'       : convert_arc,
                    'LINE'      : convert_line,
                    'ELLIPSE'   : convert_ellipse,
                    'HATCH'     : convert_hatch,
                    'LWPOLYLINE': convert_lwpolyline,
               }

    converters['INSERT'] = partial(convert_insert, converters=converters)#needed to allow converter parameter
    #converters['INSERT'] = lambda x:([],True)
    
    layers = dict()
    for e in dxf.modelspace():
        objects,append = converters[e.dxftype()](e)
        for index,o in enumerate(objects):
            if o["layer"] not in layers:    layers[o["layer"]] = [o]
            elif append:                    layers[o["layer"]].append(o)
            else:                           layers[o["layer"]].insert(index,o)
            
    write_gerber(layers)

if 2 <= len(sys.argv):    
    convert_dxf2gerber(sys.argv[1])
else:
    for filename in glob.glob("*.dxf"):
        convert_dxf2gerber(filename)

