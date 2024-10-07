module SteelDeck

using CUFSM, SectionProperties, AISIS100, SDIComposite, OrderedCollections, NonlinearSolve

include("BareProperties.jl")

include("CompositeProperties.jl")

include("ConstructionSpans.jl")

end # module SteelDeck
