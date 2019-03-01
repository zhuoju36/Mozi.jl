module Enums
export SectionType,MaterialType,LoadCaseType
@enum SectionType begin
    GENERAL_SECTION=0
    ISECTION=1
    HSECTION=2
    BOX=3
    PIPE=4
    CIRCLE=5
    RECTANGLE=6
end

@enum LoadCaseType begin
    STATIC
    MODAL
    BUCKLING
    TIME_HISTORY
    SPECTRUM
end

@enum MaterialType begin
    GENERAL_MATERIAL
    ISOELASTIC
    UNIAXIAL_METAL
end

end
