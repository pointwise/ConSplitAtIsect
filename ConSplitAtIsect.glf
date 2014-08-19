#
# Copyright 2012 (c) Pointwise, Inc.
# All rights reserved.
#
# This sample script is not supported by Pointwise, Inc.
# It is provided freely for demonstration purposes only.
# SEE THE WARRANTY DISCLAIMER AT THE BOTTOM OF THIS FILE.
#

package require PWI_Glyph 2

# enable Tk in Glyph 2
pw::Script loadTk

set isectTol 1e-5

set conName(1) ""
set conName(2) ""

#
# If there are exactly two connectors selected, initialize the Tk widget with
# the selection. Otherwise, leave the Tk widget empty.
#
set mask [pw::Display createSelectionMask -requireConnector {Dimensioned}]
pw::Display getSelectedEntities -selectionmask $mask selection

if {[llength $selection(Connectors)] == 2} {
  set conName(1) [[lindex $selection(Connectors) 0] getName]
  set conName(2) [[lindex $selection(Connectors) 1] getName]
}

proc getConPtData { con arc } {
  return [list $arc [$con getXYZ -arc $arc]]
}

proc getNextDeltaArc { con sLow sHigh percentTol ptTol } {
  set pLow [$con getXYZ -arc $sLow]
  set doit 1
  set s $sHigh
  set sBase [$con getLength -arc 1.0]
  set tol [expr {1.0 + $percentTol}]

  while { $doit } {
    set doit 0
    set p [$con getXYZ -arc $s]

    set pDelta [pwu::Vector3 length [pwu::Vector3 subtract $p $pLow]]

    set sDelta [expr {$sBase * ($s - $sLow)}]
    if {$sDelta >= $tol * $pDelta} {
      set doit 1
    } else {
      set pMid [pwu::Vector3 add $pLow \
        [pwu::Vector3 scale [pwu::Vector3 subtract $p $pLow] 0.5]]

      set sMid [expr {$sLow + 0.5 * ($s - $sLow)}]
      set vMid [$con getXYZ -arc $sMid]

      set vDelta [pwu::Vector3 length [pwu::Vector3 subtract $vMid $pMid]]

      if {$vDelta > $ptTol} {
        set doit 1
      }
    }

    if {$doit} {
      set s [expr {$sLow + 0.75 * ($s - $sLow)}]
    }
  }
  return $s
}

proc getConApprox { con sLow sHigh percentTol ptTol } {
  global isectTol

  set conData [list [getConPtData $con $sLow]]
  set sCurr $sLow
  while {$sCurr < $sHigh} {
    set sNext [getNextDeltaArc $con $sCurr $sHigh $percentTol $ptTol]
    lappend conData [getConPtData $con $sNext]
    set sCurr $sNext
  }
  return $conData
}

#
# Given 2 segments defined by their endpoints so that
#
#   p = p0 + s * (p1 - p0)  for 0 <= s <= 1
#   q = q1 + t * (q1 - q0)  for 0 <= t <= 1
#
# the minimum distance occurs when the line between p and q is perpendicular
# to both segments.
#
#   (p - q) . (p1 - p0) = 0
#   (p - q) . (q1 - q0) = 0
#
proc minDistBetweenSegs { p0 p1 q0 q1 } {

  set pd [pwu::Vector3 subtract $p1 $p0]
  set qd [pwu::Vector3 subtract $q1 $q0]
  set pq [pwu::Vector3 subtract $p0 $q0]
  set A [pwu::Vector3 dot $pd $pd]
  set B [pwu::Vector3 dot $pd $qd]
  set C [pwu::Vector3 dot $pq $pd]
  set E [pwu::Vector3 dot $qd $qd]
  set F [pwu::Vector3 dot $pq $qd]

  # Den should always be non-negative:
  #   |pd|^2 * |qd|^2 - |pd||qd|cos(a) * |pd||qd|cos(a) >= 0
  set den [expr {$A * $E - $B * $B}]

  # If it's close to zero, the lines are parallel and the distance is constant
  if {$den < 1.e-08} {
    set s 0.0
    if {$E < 1.e-08} {
      set t 0.0
    } else {
      set t [expr {$F / $E}]
    }
  } else {
    set s [expr {double($B * $F - $C * $E) / $den}]
    set t [expr {double($A * $F - $C * $B) / $den}]
  }
  # Clamp to the segments
  if {$s < 0.0 } {
    set p $p0
  } elseif {$s > 1.0} {
    set p $p1
  } else {
    set p [pwu::Vector3 add $p0 [pwu::Vector3 scale $pd $s]]
  }
  if {$t < 0.0 } {
    set q $q0
  } elseif {$t > 1.0} {
    set q $q1
  } else {
    set q [pwu::Vector3 add $q0 [pwu::Vector3 scale $qd $t]]
  }
  # Calculate the squared distance between the points
  set delta [pwu::Vector3 subtract $p $q]

  return [list $p $q [pwu::Vector3 dot $delta $delta] $s $t]
}

proc getClosestPtFromApprox { cn1Data cn2Data } {
  global isectTol
  set tol [expr {1.1 * [pw::Grid getNodeTolerance]}]

  set tol2 [expr {$tol * $tol}]
  set iMax [llength $cn1Data]
  set jMax [llength $cn2Data]
  set s1 [lindex [lindex $cn1Data 0] 0]
  set p1 [lindex [lindex $cn1Data 0] 1]
  set minDist $isectTol
  set iMin -1
  set jMin -1
  set pMin $p1
  set qMin [lindex [lindex $cn2Data 0] 1]
  for {set i 1} {$i < $iMax} {incr i} {
    set s0 $s1
    set p0 $p1
    set s1 [lindex [lindex $cn1Data $i] 0]
    set p1 [lindex [lindex $cn1Data $i] 1]
    set t1 [lindex [lindex $cn2Data 0] 0]
    set q1 [lindex [lindex $cn2Data 0] 1]
    for {set j 1} {$j < $jMax} {incr j} {
      set t0 $t1
      set q0 $q1
      set t1 [lindex [lindex $cn2Data $j] 0]
      set q1 [lindex [lindex $cn2Data $j] 1]
      set result [minDistBetweenSegs $p0 $p1 $q0 $q1]
      set dist [lindex $result 2]
      if {$dist < $minDist} {
        set useIt 1
        set p [lindex $result 0]
        set q [lindex $result 1]
        set s [lindex $result 3]
        set t [lindex $result 4]
        # Check to see if intersection is too close to an end
        if {0.0 >= $s0} {
          set d [pwu::Vector3 subtract $p $p0]

          if { [pwu::Vector3 dot $d $d] < $tol2 } {
            set useIt 0
          }
        }
        if {1.0 <= $s1} {
          set d [pwu::Vector3 subtract $p $p1]

          if { [pwu::Vector3 dot $d $d] < $tol2} {
            set useIt 0
          }
        }
        if {0.0 >= $t0} {
          set d [pwu::Vector3 subtract $q $q0]

          if { [pwu::Vector3 dot $d $d] < $tol2} {
            set useIt 0
          }
        }
        if {1.0 <= $t1} {
          set d [pwu::Vector3 subtract $q $q1]

          if { [pwu::Vector3 dot $d $d] < $tol2} {
            set useIt 0
          }
        }
        if {$useIt} {
          set minDist $dist
          set iMin [expr {$i - 1}]
          set jMin [expr {$j - 1}]
          set pMin $p
          set qMin $q
        }
      }
    }
  }

  return [list [expr {sqrt($minDist)}] $iMin $jMin $pMin $qMin]
}

proc findIntersection { cn1 cn2 } {
  global isectTol
  set solutions [list]
  set tol [expr {100.0 * $isectTol}]
  set cn1Data [getConApprox $cn1 0.0 1.0 0.001 $tol]
  set cn2Data [getConApprox $cn2 0.0 1.0 0.001 $tol]
  set result [getClosestPtFromApprox $cn1Data $cn2Data]

  set minDist [lindex $result 0]
  set iMin [lindex $result 1]
  set jMin [lindex $result 2]
  set pMin [lindex $result 3]
  set qMin [lindex $result 4]
  if {$iMin > -1} {
    set iMax [llength $cn1Data]
    if {$iMin + 2 < $iMax} {
      set iMax [expr {$iMin + 2}]
    } else {
      set iMax [expr {$iMax - 1}]
    }
    if {$iMin > 0} {
      set iMin [expr {$iMin - 1}]
    }
    set jMax [llength $cn2Data]
    if {$jMin + 2 < $jMax} {
      set jMax [expr {$jMin + 2}]
    } else {
      set jMax [expr {$jMax - 1}]
    }
    if {$jMin > 0} {
      set jMin [expr {$jMin - 1}]
    }
    set s1Min [lindex [lindex $cn1Data $iMin] 0]
    set s1Max [lindex [lindex $cn1Data $iMax] 0]
    set s2Min [lindex [lindex $cn2Data $jMin] 0]
    set s2Max [lindex [lindex $cn2Data $jMax] 0]
    set tol $isectTol
    set cn1Data [getConApprox $cn1 $s1Min $s1Max 0.0001 $tol]
    set cn2Data [getConApprox $cn2 $s2Min $s2Max 0.0001 $tol]
    set result [getClosestPtFromApprox $cn1Data $cn2Data]
    set minDist [lindex $result 0]
    set iMin [lindex $result 1]
    set jMin [lindex $result 2]
    if {$iMin > -1} {
      if {$minDist <= $isectTol && $minDist <= [pw::Grid getNodeTolerance]} {
        set sol1 [lindex $result 3]
        set sol2 [lindex $result 4]
        set solutions [list \
          [pwu::Vector3 add $sol1 [pwu::Vector3 scale \
            [pwu::Vector3 subtract $sol2 $sol1] 0.5]]]

      } else {
        set solutions [list [lindex $result 3] [lindex $result 4]]
      }
    }
  }
  return $solutions
}

proc validTolerance { tol } {
  return [expr [string is double $tol] && $tol > 0.0]
}

proc splitConnectors { } {
  global conName isectTol

  set nodeTol [pw::Grid getNodeTolerance]

  foreach i { 1 2 } {
    if [catch { pw::Grid getByName $conName($i) } con($i)] {
      tk_messageBox -icon error -message \
          "$conName($i) is not a valid connector name" -parent . \
          -title "Split Error" -type ok
      return
    }
  }

  set inters [findIntersection $con(1) $con(2)]

  if { [llength $inters] == 0 } {
    tk_messageBox -icon error -message \
        "Connectors do not intersect within the given tolerance." -parent . \
        -title "No Intersection Found" -type ok
    return
  }

  set inter(1) [lindex $inters 0]
  if {1 < [llength $inters]} {
    set inter(2) [lindex $inters 1]
  } else {
    set inter(2) $inter(1)
  }

  set errors ""

  foreach i { 1 2 } {
    if { [pwu::Vector3 length [pwu::Vector3 subtract $inter($i) \
         [$con($i) getXYZ -arc 0]]] <= $nodeTol || \
         [pwu::Vector3 length [pwu::Vector3 subtract $inter($i) \
         [$con($i) getXYZ -arc 1]]] <= $nodeTol} {

      set errors "${errors}Skipping intersection on $conName(1) because it is \
                  too close to an end point\n"
    } else {
      set pt $inter($i)
        if [catch {
            set tc [$con($i) split [$con($i) getParameter -closest $pt]]
            pw::Display update
          } msg] {
            if { ! [string equal $msg \
                "ERROR: point is not on the connector\n"]} {
              set errors "$errors$msg"
            }
          }
    }
  }

  if [string length $errors] {
    tk_messageBox -icon error -message \
        "Errors occured during split command:\n$errors" -parent . \
        -title "Split Error" -type ok
  } else {
    puts "$conName(1) and $conName(2) were split at the intersection."
  }
}

######################################################################
#  PROC: CreateLabelFrame
#    Creates a fancy label frame widget
#    Returns the new frame
#
proc CreateLabelFrame {w args} {
  #-- strip extraneous '.'s in window name
  set w [string trim $w "."]
  set w ".$w"
  frame $w -bd 0
  label $w.l
  frame $w.f -bd 2 -relief groove
  frame $w.f.spc -height 5
  pack $w.f.spc
  frame $w.f.f
  pack $w.f.f
  set text {}
  set font {}
  set padx 3
  set pady 7
  set ipadx 2
  set ipady 9
  set ipady 5
  foreach {tag value} $args {
    switch -- $tag {
      -font  {set font $value}
      -text  {set text $value}
      -padx  {set padx $value}
      -pady  {set pady $value}
      -ipadx {set ipadx $value}
      -ipady {set ipady $value}
      -bd     {$w.f config -bd $value}
      -relief {$w.f config -relief $value}
    }
  }
  if {"$font"!=""} {
    $w.l config -font $font
  }
  $w.l config -text $text
  pack $w.f -padx $padx -pady $pady -fill both -expand 1
  place $w.l -x [expr $padx+10] -y $pady -anchor w
  pack $w.f.f -padx $ipadx -pady $ipady -fill both -expand 1
  raise $w.l
  return $w.f.f
}

proc doIt { } {
  global isectTol conName

  if [checkDoIt] {
    set oldCursor [. cget -cursor]
    . configure -cursor watch
    update
    splitConnectors
    . configure -cursor $oldCursor
    set conName(1) ""
    set conName(2) ""
    checkDoIt
  } else {
    tk_messageBox -icon error -message \
        "Invalid tolerance ($isectTol), could not split connectors." \
        -parent . -title "Tolerance Error" -type ok
  }
}

proc selectConnectors { } {
  global conName

  if { [llength [pw::Grid getAll -type pw::Connector]] > 1 } {
    wm withdraw .

    set mask [pw::Display createSelectionMask -requireConnector {Dimensioned}]
    set exclude [list]

    foreach i { 1 2 } {
      if { $i == 1 } {
        pw::Display selectEntities -selectionmask $mask \
            -description "Select first connector" -single \
            selection
      } else {
        pw::Display selectEntities -selectionmask $mask \
            -description "Select second connector" -single -exclude $exclude \
            selection
      }

      if [llength $selection(Connectors)] {
        set con [lindex $selection(Connectors) 0]
        set conName($i) [$con getName]
        lappend exclude $con
      } else {
        set conName($i) ""
        break
      }
    }

    if { ! [string length $conName(1)] || ! [string length $conName(2)] } {
      set conName(1) ""
      set conName(2) ""
    }

    if {[winfo exists .]} {
      wm deiconify .
    }
  } else {
    tk_messageBox -icon error \
      -message "No connectors are available for picking" -parent . \
      -title "No Valid Entities" -type ok
  }

  checkDoIt
}

proc validateConName { name field } {
  if [catch { pw::Grid getByName $name }] {
    $field configure -background #FFCCCC
  } else {
    $field configure -background #FFFFFF
  }
  after idle checkDoIt
  return 1
}

proc validateTolerance { tol field } {
  if [validTolerance $tol] {
    $field configure -background #FFFFFF
  } else {
    $field configure -background #FFCCCC
  }
  after idle checkDoIt
  return 1
}

proc checkDoIt { } {
  global conName isectTol

  set canDoIt 1

  if { ! [validTolerance $isectTol] } {
    set canDoIt 0
  } else {
    foreach i { 1 2 } {
      if [catch { pw::Grid getByName $conName($i) } con($i)] {
        set canDoIt 0
        break
      }
    }
  }

  if $canDoIt {
    .buttons.apply configure -state normal
    .buttons.ok configure -state normal
  } else {
    .buttons.apply configure -state disabled
    .buttons.ok configure -state disabled
  }

  return $canDoIt
}

proc makeWindow { } {
  global conName isectTol

  label .title -text "Split Two Connectors at Intersection"
  set font [.title cget -font]
  set bold [font create -family [font actual $font -family] -weight bold]

  .title configure -font $bold
  pack .title -expand 1 -side top

  pack [frame .hr1 -bd 1 -height 2 -relief sunken] -fill x -pady 2
  pack [frame .content] -fill both -side top -padx 2

  set frmCons [CreateLabelFrame .connectors -text "Connectors"]
  pack [button .pick -command { selectConnectors } -text "Pick" -width 5] \
     -side right -padx 5 -in $frmCons

  pack [frame .con(1)] -in $frmCons -side top -pady 2
  pack [label .con(1).lbl -text "#1:"] -side left -padx 3
  pack [entry .con(1).name -textvariable conName(1) \
     -validatecommand {validateConName %P %W} -validate key] \
     -side right -padx 7

  pack [frame .con(2)] -in $frmCons -side bottom -pady 2
  pack [label .con(2).lbl -text "#2:"] -side left -padx 3
  pack [entry .con(2).name -textvariable conName(2) \
     -validatecommand {validateConName %P %W} -validate key] \
     -side right -padx 7

  pack .connectors -side top -pady 3 -padx 5 -fill x

  pack [frame .tol] -side top -fill x
  pack [label .tol.lbl -text "Tolerance:"] -side left -padx 8
  pack [entry .tol.name -width 10 -textvariable isectTol] \
     -side left -padx 0

  pack [frame .hr2 -bd 1 -height 2 -relief sunken] -fill x -pady 2
  pack [frame .buttons] -fill both -side top

  pack [button .buttons.apply -text "Apply" -command { doIt } -width 5] \
    -padx 6 -side right
  pack [button .buttons.cancel -text "Close" -command { exit } -width 5] \
    -padx 3 -side right
  pack [button .buttons.ok -text "OK" -command { doIt; exit } -width 5] \
    -padx 3 -side right

  pack [label .buttons.logo -image [pwLogo] -bd 0 -relief flat] \
     -side left -padx 5

  .tol.name configure -validatecommand {validateTolerance %P %W} -validate key

  wm title . "Split Connectors at Intersection"

  checkDoIt
}

proc pwLogo { } {
  set logoData "
R0lGODlheAAYAIcAAAAAAAICAgUFBQkJCQwMDBERERUVFRkZGRwcHCEhISYmJisrKy0tLTIyMjQ0
NDk5OT09PUFBQUVFRUpKSk1NTVFRUVRUVFpaWlxcXGBgYGVlZWlpaW1tbXFxcXR0dHp6en5+fgBi
qQNkqQVkqQdnrApmpgpnqgpprA5prBFrrRNtrhZvsBhwrxdxsBlxsSJ2syJ3tCR2siZ5tSh6tix8
ti5+uTF+ujCAuDODvjaDvDuGujiFvT6Fuj2HvTyIvkGKvkWJu0yUv2mQrEOKwEWNwkaPxEiNwUqR
xk6Sw06SxU6Uxk+RyVKTxlCUwFKVxVWUwlWWxlKXyFOVzFWWyFaYyFmYx16bwlmZyVicyF2ayFyb
zF2cyV2cz2GaxGSex2GdymGezGOgzGSgyGWgzmihzWmkz22iymyizGmj0Gqk0m2l0HWqz3asznqn
ynuszXKp0XKq1nWp0Xaq1Hes0Xat1Hmt1Xyt0Huw1Xux2IGBgYWFhYqKio6Ojo6Xn5CQkJWVlZiY
mJycnKCgoKCioqKioqSkpKampqmpqaurq62trbGxsbKysrW1tbi4uLq6ur29vYCu0YixzYOw14G0
1oaz14e114K124O03YWz2Ie12oW13Im10o621Ii22oi23Iy32oq52Y252Y+73ZS51Ze81JC625G7
3JG825K83Je72pW93Zq92Zi/35G+4aC90qG+15bA3ZnA3Z7A2pjA4Z/E4qLA2KDF3qTA2qTE3avF
36zG3rLM3aPF4qfJ5KzJ4LPL5LLM5LTO4rbN5bLR6LTR6LXQ6r3T5L3V6cLCwsTExMbGxsvLy8/P
z9HR0dXV1dbW1tjY2Nra2tzc3N7e3sDW5sHV6cTY6MnZ79De7dTg6dTh69Xi7dbj7tni793m7tXj
8Nbk9tjl9N3m9N/p9eHh4eTk5Obm5ujo6Orq6u3t7e7u7uDp8efs8uXs+Ozv8+3z9vDw8PLy8vL0
9/b29vb5+/f6+/j4+Pn6+/r6+vr6/Pn8/fr8/Pv9/vz8/P7+/gAAACH5BAMAAP8ALAAAAAB4ABgA
AAj/AP8JHEiwoMGDCBMqXMiwocOHECNKnEixosWLGDNqZCioo0dC0Q7Sy2btlitisrjpK4io4yF/
yjzKRIZPIDSZOAUVmubxGUF88Aj2K+TxnKKOhfoJdOSxXEF1OXHCi5fnTx5oBgFo3QogwAalAv1V
yyUqFCtVZ2DZceOOIAKtB/pp4Mo1waN/gOjSJXBugFYJBBflIYhsq4F5DLQSmCcwwVZlBZvppQtt
D6M8gUBknQxA879+kXixwtauXbhheFph6dSmnsC3AOLO5TygWV7OAAj8u6A1QEiBEg4PnA2gw7/E
uRn3M7C1WWTcWqHlScahkJ7NkwnE80dqFiVw/Pz5/xMn7MsZLzUsvXoNVy50C7c56y6s1YPNAAAC
CYxXoLdP5IsJtMBWjDwHHTSJ/AENIHsYJMCDD+K31SPymEFLKNeM880xxXxCxhxoUKFJDNv8A5ts
W0EowFYFBFLAizDGmMA//iAnXAdaLaCUIVtFIBCAjP2Do1YNBCnQMwgkqeSSCEjzzyJ/BFJTQfNU
WSU6/Wk1yChjlJKJLcfEgsoaY0ARigxjgKEFJPec6J5WzFQJDwS9xdPQH1sR4k8DWzXijwRbHfKj
YkFO45dWFoCVUTqMMgrNoQD08ckPsaixBRxPKFEDEbEMAYYTSGQRxzpuEueTQBlshc5A6pjj6pQD
wf9DgFYP+MPHVhKQs2Js9gya3EB7cMWBPwL1A8+xyCYLD7EKQSfEF1uMEcsXTiThQhmszBCGC7G0
QAUT1JS61an/pKrVqsBttYxBxDGjzqxd8abVBwMBOZA/xHUmUDQB9OvvvwGYsxBuCNRSxidOwFCH
J5dMgcYJUKjQCwlahDHEL+JqRa65AKD7D6BarVsQM1tpgK9eAjjpa4D3esBVgdFAB4DAzXImiDY5
vCFHESko4cMKSJwAxhgzFLFDHEUYkzEAG6s6EMgAiFzQA4rBIxldExBkr1AcJzBPzNDRnFCKBpTd
gCD/cKKKDFuYQoQVNhhBBSY9TBHCFVW4UMkuSzf/fe7T6h4kyFZ/+BMBXYpoTahB8yiwlSFgdzXA
5JQPIDZCW1FgkDVxgGKCFCywEUQaKNitRA5UXHGFHN30PRDHHkMtNUHzMAcAA/4gwhUCsB63uEF+
bMVB5BVMtFXWBfljBhhgbCFCEyI4EcIRL4ChRgh36LBJPq6j6nS6ISPkslY0wQbAYIr/ahCeWg2f
ufFaIV8QNpeMMAkVlSyRiRNb0DFCFlu4wSlWYaL2mOp13/tY4A7CL63cRQ9aEYBT0seyfsQjHedg
xAG24ofITaBRIGTW2OJ3EH7o4gtfCIETRBAFEYRgC06YAw3CkIqVdK9cCZRdQgCVAKWYwy/FK4i9
3TYQIboE4BmR6wrABBCUmgFAfgXZRxfs4ARPPCEOZJjCHVxABFAA4R3sic2bmIbAv4EvaglJBACu
IxAMAKARBrFXvrhiAX8kEWVNHOETE+IPbzyBCD8oQRZwwIVOyAAXrgkjijRWxo4BLnwIwUcCJvgP
ZShAUfVa3Bz/EpQ70oWJC2mAKDmwEHYAIxhikAQPeOCLdRTEAhGIQKL0IMoGTGMgIBClA9QxkA3U
0hkKgcy9HHEQDcRyAr0ChAWWucwNMIJZ5KilNGvpADtt5JrYzKY2t8nNbnrzm+B8SEAAADs="

  return [image create photo -format GIF -data $logoData]
}

makeWindow
::tk::PlaceWindow . widget
tkwait window .

#
# DISCLAIMER:
# TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW, POINTWISE DISCLAIMS
# ALL WARRANTIES, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED
# TO, IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE, WITH REGARD TO THIS SCRIPT.  TO THE MAXIMUM EXTENT PERMITTED
# BY APPLICABLE LAW, IN NO EVENT SHALL POINTWISE BE LIABLE TO ANY PARTY
# FOR ANY SPECIAL, INCIDENTAL, INDIRECT, OR CONSEQUENTIAL DAMAGES
# WHATSOEVER (INCLUDING, WITHOUT LIMITATION, DAMAGES FOR LOSS OF
# BUSINESS INFORMATION, OR ANY OTHER PECUNIARY LOSS) ARISING OUT OF THE
# USE OF OR INABILITY TO USE THIS SCRIPT EVEN IF POINTWISE HAS BEEN
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGES AND REGARDLESS OF THE
# FAULT OR NEGLIGENCE OF POINTWISE.
#
