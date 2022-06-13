import cobra

def get_reactions_to_convert(model):
    reactions_not_to_convert = []
    all_reactions = model.reactions
    reactions_to_convert = []

    transport_reactions = []

    exchanges = model.exchanges

    demands = model.demands
    sinks = model.sinks

    exchange_reactions = []
    demand_reactions = []
    sink_reactions = []
    exchange_reactions.extend(exchanges)
    demand_reactions.extend(demands)
    sink_reactions.extend(sink_reactions)

    for group in model.groups:
        reactions_list = []
        if "transport" in str(group.name).lower() or "drain" in str(group.name).lower():
            for member in group.members:
                reactions_list.append(member)
            transport_reactions.extend(reactions_list)

    for reaction in all_reactions:
        reaction_id = reaction.id

        add = True

        # exclude exchange reactions from the comparison
        if len(exchanges) > 0 and add:
            if reaction in exchanges:
                add = False

        if "drain" in str(reaction.name).lower() and add:
            exchange_reactions.append(reaction)
            print(reaction_id)

        if "exchange" in str(reaction.name).lower() and add:
            exchange_reactions.append(reaction)

        if str(reaction_id).upper().startswith("EX_") and add:
            exchange_reactions.append(reaction)

        # exclude demand reactions from the comparison
        if len(demands) > 0 and add:
            if reaction in demands:
                add = False

        # exclude transport reactions from the comparison
        if len(sinks) > 0 and add:
            if reaction in sinks:
                pass

        # exclude transport reactions from the comparison
        if "transport" in str(reaction.name).lower() and add:
            transport_reactions.append(reaction)

        if str(reaction_id).upper().startswith("TRANS-RXN") and add:
            transport_reactions.append(reaction)

        if "TRANS-RXN" in str(reaction.name).upper() and add:
            transport_reactions.append(reaction)


    print("Transport reactions: %d" % len(transport_reactions))
    print("Sink reactions: %d" % len(sink_reactions))
    print("demand reactions: %d" % len(demand_reactions))
    print("exchange reactions: %d" % len(exchange_reactions))
    print("compartments: %d" % len(model.compartments))

if __name__ == "__main__":
    model = cobra.io.sbml.read_sbml_model("C:/Users/maria/Documents/GitHub/Project_DraftModels/DraftModels/Models/AuReMe/models/aureme_model.sbml")
    get_reactions_to_convert(model)
    model = cobra.io.sbml.read_sbml_model("C:/Users/maria/Documents/GitHub/Project_DraftModels/DraftModels/Models/CarveMe/models/model_carveme.xml")
    get_reactions_to_convert(model)
    model = cobra.io.sbml.read_sbml_model("C:/Users/maria/Documents/GitHub/Project_DraftModels/DraftModels/Models/KBase/model/kbase_model.xml")
    get_reactions_to_convert(model)
    #model = cobra.io.sbml.read_sbml_model("C:/Users/maria/Documents/GitHub/Project_DraftModels/DraftModels/Models/")
    #get_reactions_to_convert(model)
    model = cobra.io.sbml.read_sbml_model("C:/Users/maria/Documents/GitHub/Project_DraftModels/DraftModels/Models/merlin(BIT)/model/BIT_model.xml")
    get_reactions_to_convert(model)