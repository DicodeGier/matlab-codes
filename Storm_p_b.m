Function ErlangB (E As Double, m As Integer) As Double
    Dim InvB As Double
    Dim j As Integer

    InvB = 1.0
    For j = 1 To m
        InvB = 1.0 + InvB * j / E 
    Next j
    ErlangB = 1.0 / InvB
End Function