LIST ( APPEND SConfig_files
CMakeLists.txt
ConfigException.hh
ConfigFileReader.cxx
ConfigFileReader.hh
ConfigMap.hh
ConfigObject.cxx
ConfigObject.hh
ConfigOption.hh
Counter.hh
Factory.hh
GenericOption.hh
MMap.hh
NamedObject.hh
ObjectProvider.hh
OptionValidation.hh
Provider.hh
SConfigPtr.hh
SharedPtr.hh
StringManip.hh
)

SF_ADD_KERNEL_LIBRARY ( SConfig )

