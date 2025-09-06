#include <bayeux/dpp/chain_module.h>

#include <falaise/property_set.h>

#include <cstdio>
#include <limits>

////////////////////////////////////////////////////////////////

class bank_removal_module : public dpp::chain_module
{
public:
  bank_removal_module();
  ~bank_removal_module();

  void initialize(const datatools::properties &,
		  datatools::service_manager &,
		  dpp::module_handle_dict_type &);
  
  dpp::chain_module::process_status process(datatools::things &);

  void finalize();

private:
  std::vector<std::string> _banks_to_remove_;

  DPP_MODULE_REGISTRATION_INTERFACE(bank_removal_module);
};

DPP_MODULE_REGISTRATION_IMPLEMENT(bank_removal_module, "bank_removal_module");

////////////////////////////////////////////////////////////////

bank_removal_module::bank_removal_module() : dpp::chain_module()
{
  
}

bank_removal_module::~bank_removal_module()
{
  this->finalize();
}

void bank_removal_module::initialize(const datatools::properties & config,
			   datatools::service_manager & ,
			   dpp::module_handle_dict_type & )
{
  const falaise::property_set fps {config};

  if (config.has_key("banks_to_remove"))
    config.fetch("banks_to_remove", _banks_to_remove_);

  this->_set_initialized(true);
}

dpp::chain_module::process_status bank_removal_module::process(datatools::things &event)
{
  for (const std::string & bank_to_remove : _banks_to_remove_)
    // if (event.has(bank_to_remove))
      event.remove(bank_to_remove);

  return dpp::base_module::PROCESS_SUCCESS;
}

void bank_removal_module::finalize()
{
  this->_set_initialized(false);
}
