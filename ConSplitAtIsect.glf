#############################################################################
#
# (C) 2021 Cadence Design Systems, Inc. All rights reserved worldwide.
#
# This sample script is not supported by Cadence Design Systems, Inc.
# It is provided freely for demonstration purposes only.
# SEE THE WARRANTY DISCLAIMER AT THE BOTTOM OF THIS FILE.
#
#############################################################################

package require PWI_Glyph 2

# enable Tk in Glyph 2
pw::Script loadTk

set isectTol 1e-5

set conName(1) ""
set conName(2) ""

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

  pack [label .buttons.logo -image [cadenceLogo] -bd 0 -relief flat] \
     -side left -padx 5

  .tol.name configure -validatecommand {validateTolerance %P %W} -validate key

  wm title . "Split Connectors at Intersection"

  checkDoIt
}

proc cadenceLogo { } {
  set logoData "
R0lGODlhgAAYAPQfAI6MjDEtLlFOT8jHx7e2tv39/RYSE/Pz8+Tj46qoqHl3d+vq62ZjY/n4+NT
T0+gXJ/BhbN3d3fzk5vrJzR4aG3Fubz88PVxZWp2cnIOBgiIeH769vtjX2MLBwSMfIP///yH5BA
EAAB8AIf8LeG1wIGRhdGF4bXD/P3hwYWNrZXQgYmVnaW49Iu+7vyIgaWQ9Ilc1TTBNcENlaGlIe
nJlU3pOVGN6a2M5ZCI/PiA8eDp4bXBtdGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1w
dGs9IkFkb2JlIFhNUCBDb3JlIDUuMC1jMDYxIDY0LjE0MDk0OSwgMjAxMC8xMi8wNy0xMDo1Nzo
wMSAgICAgICAgIj48cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudy5vcmcvMTk5OS8wMi
8yMi1yZGYtc3ludGF4LW5zIyI+IDxyZGY6RGVzY3JpcHRpb24gcmY6YWJvdXQ9IiIg/3htbG5zO
nhtcE1NPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvbW0vIiB4bWxuczpzdFJlZj0iaHR0
cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL3NUcGUvUmVzb3VyY2VSZWYjIiB4bWxuczp4bXA9Imh
0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtcE1NOk9yaWdpbmFsRG9jdW1lbnRJRD0idX
VpZDoxMEJEMkEwOThFODExMUREQTBBQzhBN0JCMEIxNUM4NyB4bXBNTTpEb2N1bWVudElEPSJ4b
XAuZGlkOkIxQjg3MzdFOEI4MTFFQjhEMv81ODVDQTZCRURDQzZBIiB4bXBNTTpJbnN0YW5jZUlE
PSJ4bXAuaWQ6QjFCODczNkZFOEI4MTFFQjhEMjU4NUNBNkJFRENDNkEiIHhtcDpDcmVhdG9yVG9
vbD0iQWRvYmUgSWxsdXN0cmF0b3IgQ0MgMjMuMSAoTWFjaW50b3NoKSI+IDx4bXBNTTpEZXJpZW
RGcm9tIHN0UmVmOmluc3RhbmNlSUQ9InhtcC5paWQ6MGE1NjBhMzgtOTJiMi00MjdmLWE4ZmQtM
jQ0NjMzNmNjMWI0IiBzdFJlZjpkb2N1bWVudElEPSJ4bXAuZGlkOjBhNTYwYTM4LTkyYjItNDL/
N2YtYThkLTI0NDYzMzZjYzFiNCIvPiA8L3JkZjpEZXNjcmlwdGlvbj4gPC9yZGY6UkRGPiA8L3g
6eG1wbWV0YT4gPD94cGFja2V0IGVuZD0iciI/PgH//v38+/r5+Pf29fTz8vHw7+7t7Ovp6Ofm5e
Tj4uHg397d3Nva2djX1tXU09LR0M/OzczLysnIx8bFxMPCwcC/vr28u7q5uLe2tbSzsrGwr66tr
KuqqainpqWko6KhoJ+enZybmpmYl5aVlJOSkZCPjo2Mi4qJiIeGhYSDgoGAf359fHt6eXh3dnV0
c3JxcG9ubWxramloZ2ZlZGNiYWBfXl1cW1pZWFdWVlVUU1JRUE9OTUxLSklIR0ZFRENCQUA/Pj0
8Ozo5ODc2NTQzMjEwLy4tLCsqKSgnJiUkIyIhIB8eHRwbGhkYFxYVFBMSERAPDg0MCwoJCAcGBQ
QDAgEAACwAAAAAgAAYAAAF/uAnjmQpTk+qqpLpvnAsz3RdFgOQHPa5/q1a4UAs9I7IZCmCISQwx
wlkSqUGaRsDxbBQer+zhKPSIYCVWQ33zG4PMINc+5j1rOf4ZCHRwSDyNXV3gIQ0BYcmBQ0NRjBD
CwuMhgcIPB0Gdl0xigcNMoegoT2KkpsNB40yDQkWGhoUES57Fga1FAyajhm1Bk2Ygy4RF1seCjw
vAwYBy8wBxjOzHq8OMA4CWwEAqS4LAVoUWwMul7wUah7HsheYrxQBHpkwWeAGagGeLg717eDE6S
4HaPUzYMYFBi211FzYRuJAAAp2AggwIM5ElgwJElyzowAGAUwQL7iCB4wEgnoU/hRgIJnhxUlpA
SxY8ADRQMsXDSxAdHetYIlkNDMAqJngxS47GESZ6DSiwDUNHvDd0KkhQJcIEOMlGkbhJlAK/0a8
NLDhUDdX914A+AWAkaJEOg0U/ZCgXgCGHxbAS4lXxketJcbO/aCgZi4SC34dK9CKoouxFT8cBNz
Q3K2+I/RVxXfAnIE/JTDUBC1k1S/SJATl+ltSxEcKAlJV2ALFBOTMp8f9ihVjLYUKTa8Z6GBCAF
rMN8Y8zPrZYL2oIy5RHrHr1qlOsw0AePwrsj47HFysrYpcBFcF1w8Mk2ti7wUaDRgg1EISNXVwF
lKpdsEAIj9zNAFnW3e4gecCV7Ft/qKTNP0A2Et7AUIj3ysARLDBaC7MRkF+I+x3wzA08SLiTYER
KMJ3BoR3wzUUvLdJAFBtIWIttZEQIwMzfEXNB2PZJ0J1HIrgIQkFILjBkUgSwFuJdnj3i4pEIlg
eY+Bc0AGSRxLg4zsblkcYODiK0KNzUEk1JAkaCkjDbSc+maE5d20i3HY0zDbdh1vQyWNuJkjXnJ
C/HDbCQeTVwOYHKEJJwmR/wlBYi16KMMBOHTnClZpjmpAYUh0GGoyJMxya6KcBlieIj7IsqB0ji
5iwyyu8ZboigKCd2RRVAUTQyBAugToqXDVhwKpUIxzgyoaacILMc5jQEtkIHLCjwQUMkxhnx5I/
seMBta3cKSk7BghQAQMeqMmkY20amA+zHtDiEwl10dRiBcPoacJr0qjx7Ai+yTjQvk31aws92JZ
Q1070mGsSQsS1uYWiJeDrCkGy+CZvnjFEUME7VaFaQAcXCCDyyBYA3NQGIY8ssgU7vqAxjB4EwA
DEIyxggQAsjxDBzRagKtbGaBXclAMMvNNuBaiGAAA7"

  return [image create photo -format GIF -data $logoData]
}

makeWindow
::tk::PlaceWindow . widget
tkwait window .

#############################################################################
#
# This file is licensed under the Cadence Public License Version 1.0 (the
# "License"), a copy of which is found in the included file named "LICENSE",
# and is distributed "AS IS." TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE
# LAW, CADENCE DISCLAIMS ALL WARRANTIES AND IN NO EVENT SHALL BE LIABLE TO
# ANY PARTY FOR ANY DAMAGES ARISING OUT OF OR RELATING TO USE OF THIS FILE.
# Please see the License for the full text of applicable terms.
#
#############################################################################
