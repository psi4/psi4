Off[ElementData::notsubprop]
props = Delete[ElementData["Properties"], Position[ElementData["Properties"], "CommonCompoundNames"][[1,1]]]
Print[{"Atomic Number", "Name", "Property", "Value", "Units", "UnitsName", "Interval", "Note"}]
Do[
    Do[
        Print[
            {
                ielement,
                ElementData[ielement], 
                prop, 
                ElementData[ElementData[ielement], prop, "Value"], 
                ElementData[ElementData[ielement], prop, "Units"],
                ElementData[ElementData[ielement], prop, "UnitsName"],
                ElementData[ElementData[ielement], prop, "Interval"],
                ElementData[ElementData[ielement], prop, "Note"]
            }
        ], 
        {prop, props}
    ], 
    {ielement, 1, 118}
]
