
import collections
import datetime
import dateutil
import matplotlib
matplotlib.use('pgf')
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import sqlite3
import warnings

from pyne.data import atomic_mass

import plotsettings
plt.style.use('seaborn')
plt.rcParams.update(plotsettings.tex_fonts())

def atom_to_mass_frac(atom):
    """Convert input in atomic fraction to mass fraction"""
    if type(atom) is np.ndarray:
        if len(atom) != 6:
            raise ValueError("Array has to be of len 6 (U232 to U236, U238)")
        isotopes = [f"92{i}0000" for i in range(232, 239) if i != 237]
        masses = np.array([atomic_mass(iso) for iso in isotopes])

        mass = atom * masses
        mass /= mass.sum()
        
        return mass
    
    mass = {}
    normalisation = 0
    for key, value in atom.items():
        if key > 250 or key < 1:
            raise KeyError('Keys have to be the atomic masses, '
                           + 'not NucIds or something else.')
        val = value * key
        mass[key] = val
        normalisation += val

    for key in mass.keys():
        mass[key] /= normalisation
    
    return mass

def mass_to_atom_frac(mass):
    """Convert input in mass fraction to atomic fraction"""
    atom = {}
    normalisation = 0
    for key, value in mass.items():
        if key > 250 or key < 1:
            raise KeyError('Keys have to be the atomic masses, '
                           + 'not NucIds or something else.')
        val = value / key
        atom[key] = val
        normalisation += val

    for key in atom.keys():
        atom[key] /= normalisation
    
    return atom
    
class Plotter:
    """
    A helper class defining plotting styles, focussing especially onto
    labeling the axes, eg., in scientific notation or using dates.
    """
    def __init__(self,
                 t_init,
                 t_ref=datetime.datetime(1980,1,1),
                 yearBase=5):
        self.t_init = t_init
        self.t_ref = t_ref
        self.yearsStep = matplotlib.dates.YearLocator(base = yearBase)
        self.months = matplotlib.dates.MonthLocator()
        self.years_fmt = matplotlib.dates.DateFormatter("%Y")
        
        v = pd.__version__.split('.')[:2]
        self.pd_version = [int(vv) for vv in v]
        if (self.pd_version[0]==0 and self.pd_version[1] < 24):
            msg =  ('Your Pandas version is < 0.24.\n'
                  + 'In function pd_to_ts, using the non-recommended '
                  + 'pd.Series.values attribute.\nIn order to prevent '
                  + 'this, upgrade Pandas to version >= 0.24.')
            warnings.warn(msg, category=ImportWarning)

        return
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Convert pandas series dataframe containing pandas datetime to
    # timestamps in seconds relative to self.t_ref
    # Returns:  - np.ndarray of shape (pd_time.size,)
    def pd_to_ts(self, pd_time):     
        # divide by 10^9 to convert ns to s 
        pd_time -= self.t_ref
        
        if (self.pd_version[0]==0 and self.pd_version[1] < 24):
            pd_time = pd_time.astype('int64').values / 1e9
        else:
            pd_time = pd_time.astype('int64').to_numpy() / 1e9
        
        return pd_time
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Convert array of dtype datetime.datetime to timestamps in seconds
    # relative to self.t_ref.
    # Returns:  - np.ndarray of shape dt.shape
    def dt_to_ts(self, dt):
        dt -=  self.t_ref
        dt = np.array([t.total_seconds() for t in dt])
        
        return dt
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Convert array of timestamps in seconds relative to self.t_ref into
    # array of datetime.datetime of same shape.
    # Returns:  - np.ndarray
    def ts_to_dt(self, ts):
        arr = np.array([self.t_ref + datetime.timedelta(seconds=int(i))
                        for i in ts])
        return arr

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculate the dates with respect to t_init, going in steps
    # as defined in arr.
    # Returns:  - a numpy array containing datetime.datetime variables
    def convert_date(self, arr, t_init=None):
        if t_init is None:
            t_init = self.t_init    
        reldelta = dateutil.relativedelta.relativedelta
        return np.array([t_init + reldelta(months=i) for i in arr])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Plot years on x-axis in steps of yearsStep years.
    # Returns:  - void function, only updates axes and figure
    def set_xaxis(self, fig, ax, yearsStep=None):
        assert (yearsStep is None
                or isinstance(yearsStep, int)
                or (isinstance(yearsStep, str) 
                    and yearsStep.lower() == "self")), \
                "'yearsStep' has to be 'None' or 'self' or of type int."
        
        ax.xaxis.set_major_formatter(self.years_fmt)
        
        if isinstance(yearsStep, int):
            yearsStep = matplotlib.dates.YearLocator(base=yearsStep)
            ax.xaxis.set_major_locator(yearsStep)
        elif yearsStep is not None:
            ax.xaxis.set_major_locator(self.yearsStep)
            
        ax.format_xdata = matplotlib.dates.DateFormatter("%Y-%m")
        fig.autofmt_xdate()
        
        return
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Plot specified axis ticks in scientific notation.
    # Returns:  - void function, only updates axes
    def set_scientific(self, ax, which_axis="y", scilimits=(-2, 3)):
        assert isinstance(which_axis, str), \
            "'which_axis' has to be of type 'str'."
        which_axis = which_axis.lower()
        assert (which_axis == "y" 
                or which_axis == "x"
                or which_axis == "both"), \
                "'which_axis' has to be 'x', 'y' or 'both'."

        ax.ticklabel_format(axis=which_axis,
                            style="sci",
                            scilimits=scilimits,
                            useOffset=None,
                            useLocale=None,
                            useMathText=True)
        return


class Analyser:
    """
    A helper class that features functions to do queries, get all of the
    agents, their archetypes and agentIds, list inventories, etc...
    """
    def __init__(self, path, fname):
        self.path = path
        self.fname = fname
        
        self.connection = sqlite3.connect(
            os.path.join(path, fname+".sqlite")
        )
        self.cursor = self.connection.cursor()
        self.t_init = self.get_initial_time()
        self.duration = self.query("Duration", "Info")[0]

        # (key,value): (agentId, name), (name, agentId),
        #              (spec, agentId), (agentId, spec)
        (self.agents,
         self.names, 
         self.specs, 
         self.agentIds) = self.get_agents()
        
        return
    
    def __del__(self):
        self.connection.close()
        print("Connection to file {}.sqlite closed.".format(self.fname))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Do a query of the database
    # Returns:  - list with query results
    def query(self, selection, table, condition="", vals=(), 
              display=False):
        if condition:
            query = "SELECT {} FROM {} WHERE {}".format(selection,
                                                          table,
                                                          condition
                                                          )
            if display:
                print(query)
            self.cursor.execute(query, vals)
        else:
            query = "SELECT {} FROM {}".format(selection, table)
            if display:
                print(query)
            self.cursor.execute(query)
        
        return self.cursor.fetchall()
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Returns:  - all or parts of the tables present in main database
    #             as a np.ndarray
    def get_tables(self, search="", display=False):
        tables = np.array(self.query("*", "SQLite_master"))
        
        if display:
            tables_ordered = tables[tables[:, 1].argsort()]
            print("\n\n\n")
            for t in tables_ordered:
                print("\n", t)
         
        return tables

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Query initial time of simulation from database
    # Returns:  - datetime.datetime with initial date
    def get_initial_time(self):
        year, month = self.query("InitialYear, InitialMonth", "Info")[0]
        return datetime.datetime(year=year, month=month, day=1)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Returns:  - dict with (key, value) = (agentId, name)
    #           - dict with (key, value) = (name, agentId)
    #           - dict with (key, value) = (spec, agentId)
    #           - dict with (key, value) = (agentId, spec)
    def get_agents(self):
        data = self.query("Spec, AgentId, Prototype", "AgentEntry")

        agents = {}
        names = {}
        specs = collections.defaultdict(list)
        agentIds = {}
        for i in range(len(data)):
            # Get the spec (without :agent: or :cycamore: etc. prefix)
            spec = str(data[i][0])
            idx = [m.start() for m in re.finditer(":", spec)][-1]
            spec = spec[idx + 1:]

            agentId = int(data[i][1])
            name = str(data[i][2])

            agents[agentId] = name
            names[name.lower()] = agentId
            specs[spec].append(agentId)
            agentIds[agentId] = spec
    
        return agents, names, specs, agentIds

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Pretty print all of the agents involved
    # Returns:  - void function
    def print_agents(self):
        i = 0
        print("AgentId   Name                     Archetype")
        print("-----------------------------------------------")
        for key, name in zip(self.agents.keys(), self.agents.values()):
            print("{:7d}   {:25}{}".format(key, name, self.agentIds[key]))
            i += 1
            if not i%5:
                print()
        return  

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Returns:  
    def get_agent_timeline(self):
        data = np.array(self.query("AgentId, Lifetime, EnterTime", 
                                   "AgentEntry"), dtype=int)
        agentId = data[:,0]
        lifetime = (np.where(data[:,1]==-1, self.duration, data[:,1])
                   - data[:,2])
        
        #DO STUFF HERE
        return agentId, lifetime
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Takes an agent name or an iterable with agent names as argument
    # Returns:  - corresponding agentId or np.ndarray with agentIds 
    def convert_to_AgentId(self, agentId):
        if isinstance(agentId, str):
            returnVal = self.names[agentId.lower()]
        
        elif isinstance(agentId, collections.abc.Iterable):
            returnVal = np.empty(len(agentId), dtype=int)
            for i, name in enumerate(agentId):
                assert isinstance(name, str), ("agentId {}".format(name)
                + "is expected to be of type 'str'")
                returnVal[i] = self.names[name.lower()]
        else:
            raise TypeError("agentId {} must be of ".format(agentId)
                            + "type 'str' or an array, list or tuple of"
                            + "datatypes 'str'")
            returnVal = -1

        return returnVal
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get net difference of transactions performed for one agent at a 
    # certain timestep. As the dictionary keys are the timesteps, each
    # timestep will be at most once in the dictionary, and it features
    # the sum of transactions performed at that timestep.
    # Returns:  - dict with (key, value) = (time, net qty) where qty > 0
    #             means incoming qty, qty < 0 means outgoing qty
    def agent_transactions(self, agentId):
        assert isinstance(agentId, (int, np.integer)), \
            "agentId has to be of type int!"
        assert agentId in self.agents, "agentId not valid!"
        
        if not (self.agentIds[agentId].lower() in ["source", "sink"]):
            msg = ("Agent '{}' with ".format(self.agents[agentId])
                   + "agentId {} is not an instance of ".format(agentId)
                   + "source or sink. Possibly wrong results because the "
                   + "destruction or conversion of material is not yet "
                   + "taken into account.")
            warnings.warn(msg)

        # quantityDict only contains the net difference at a given time
        # so the quantity has to be calculated by accumulation.
        # (key, value) = (time, net difference)
        quantity = collections.defaultdict(float)

        # Add resources created from agent.
        data = self.query("ResourceId", "ResCreators", 
                          "AgentId==?", (agentId,)) 
        res_created = tuple([i for tup in data for i in tup])
        res_created = self.query(
            "TimeCreated, Quantity",
            "Resources",
            "ResourceId in (" + ",".join(["?"]*len(res_created)) + ")",
            tuple(res_created)
        )
        for time, qty in res_created:
            quantity[time] += qty
    
        # Add resources transferred to agent.
        # (key, value) = (resId, time)
        res_received = dict(self.query(
            "ResourceId, Time",
            "Transactions",
            "ReceiverId == ?",
            (agentId,)
        ))
        if res_received:
            # (key, value) = (resId, Qty)
            transactions = dict(self.query(
                "ResourceId, Quantity",
                "Resources",
                "ResourceId in ("
                    + ",".join(["?"]*len(res_received))
                    + ")",
                tuple(res_received.keys())
            ))
            for res_id, time in res_received.items():
                quantity[time] += transactions[res_id]
    
        # Substract resources transferred from agent
        # (key, value) = (resId, time)
        res_sent = dict(self.query(
            "ResourceId, Time",
            "Transactions",
            "SenderId==?",
            (agentId,)
        ))
        if res_sent:
            # (key, value) = (resId, qty)
            transactions = dict(self.query(
                "ResourceId, Quantity",
                "Resources",
                "ResourceId in ("
                    + ",".join(["?"]*len(res_sent))
                    + ")",
                tuple(res_sent.keys())
            ))
            for res_id, time in res_sent.items():
                quantity[time] -= transactions[res_id]
    
        return quantity
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculates the inventory at each timestep for the agentId or for 
    # the list or array of agentIds that is passed to the function.
    # If fill_array is set,  missing timesteps (ie., timesteps at which 
    # no transactions have been performed) are added and the quantity 
    # value is set to the value of the previous timestep.
    # Passing an iterable can be useful if one has, eg., multiple agents
    # storing the same commodity and one wants to know the total amount
    # of said commodity.
    # Returns:  - np.array, shape(2, len(time)),
    #             where arr[0] = time, arr[1] = qty
    def agent_inventory(self, agentId, fill_array=True):
        # (key, value) = (time, qty)
        transactions = collections.defaultdict(float) 
        
        # If needed, loop over all agents.
        if isinstance(agentId, str):
            agentId = self.convert_to_AgentId(agentId)
            transactions = self.agent_transactions(agentId)
        
        elif isinstance(agentId, collections.abc.Iterable):
            if isinstance(agentId[0], str):
                agentId = self.convert_to_AgentId(agentId)
            for agent in agentId:
                trans = self.agent_transactions(agent)
                for time, net_qty in trans.items():
                    transactions[time] += net_qty
        else:
            transactions = self.agent_transactions(agentId)
        
        time_list = list(transactions.keys())
        qty_list = list(transactions.values())
        
        # Find out which timesteps are missing and append them while
        # creating an array of dim (2, max(time)+1)   
        if fill_array:
            time = np.arange(min(time_list), max(time_list)+1, 1, 
                             dtype=int)
            time_append = list(set(time) - set(time_list))
            inventory = np.array([time_list + time_append,
                                  qty_list + [0]*len(time_append)])
        else:
            inventory = np.array([time_list, qty_list])

        # Sort the array by increasing time and cumulate the commodity.
        inventory = inventory[:, inventory[0,:].argsort()]
        for i in range(1, len(inventory[1])):
            inventory[1, i] += inventory[1, i - 1]
        
        return inventory

"""
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
- - - - - - - - - - - - - - -  MAIN - - - - - - - - - - - - - - - - -
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""
if __name__=="__main__":
    import matplotlib
    matplotlib.use("TKagg")
    import matplotlib.pyplot as plt
    
    print("Starting main")  
    
    path    = "/Users/test/Uni/Masterarbeit/Pakistan/data/comprehensive"
    fname   = "300months6"
    
    test = Analyser(path, fname)
    test.print_agents()
    
    data = test.agent_inventory(27)
    
    #plt.scatter(data[0], data[1])
    #plt.show()
    




