module SteelDeck

using CUFSM, SectionProperties, AISIS100, SDIComposite, OrderedCollections, NonlinearSolve, Statistics, CrossSectionGeometry, DataFrames

include("BareDeckProperties.jl")

include("CompositeFlexuralProperties.jl")

include("WebCrippling.jl")

include("BareShear.jl")

include("OneWayCompositeShear.jl")

include("TwoWayCompositeShear.jl")

include("ConstructionLoads.jl")

include("ConstructionSpans.jl")

include("SuperimposedLoads.jl")

include("CompositeLoads.jl")

include("ConstructionLoadsUnequalSpans.jl")

include("DeckTables.jl")

end # module SteelDeck
