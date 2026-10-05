// Cutaway of the assembly from cad.py's STLs (cad.py renders it); by hand:
//   xvfb-run -a openscad -o out/cad.png --imgsize=1500,1100 --camera=0,0,34,75,0,95,300 cutaway.scad
module wedge(a0, a1) { translate([0,0,-10]) rotate([0,0,a0]) rotate_extrude(angle=a1-a0, $fn=120) square([90,100]); }
module part(name, col, a0=0, a1=0) {
  color(col) difference() { import(str("out/cad/", name, ".stl"), convexity=8); if (a1 > a0) wedge(a0, a1); }
}
// housing and stator: a 150-degree window
for (p = [["stator_core","#5b6470"], ["jacket_inner","#b9c2cc"], ["jacket_outer","#e3e7ec"], ["manifold","#d8c9a3"],
          ["drive_plate","#3f7f5f"], ["cap_rear","#aab2bb"], ["cap_front","#aab2bb"]]) part(p[0], p[1], -25, 125);
// coils: whole, except a few removed to show the bare teeth
part("coil", "#c87a2f", 60, 125);
// rotor: a narrower window, so its section and the magnets both show
for (p = [["rotor_yoke","#34404d"], ["magnet","#7aa6c2"], ["rotor_carrier","#8d99a6"]]) part(p[0], p[1], 35, 125);
part("bearing", "#e6e6e6", 35, 125); part("shaft", "#777f88");
