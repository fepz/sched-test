import os
import math
import sys
import pandas as pd
import ntpath
import re
import pickle
import matplotlib.pyplot as plt
from argparse import ArgumentParser


def test_rts(xml_file, max_fu):
    import xml.etree.cElementTree as et

    def get_rts_from_element(elem):
        """ Get the task-set from elem """

        def get_int(string):
            """ string to in int or float if it fails """
            try:
                return int(string)
            except ValueError:
                return int(float(string))

        rts_id, rts = 0, []
        if elem.tag == 'S':
            rts_id = get_int(elem.get("count"))
            for t in elem.iter("i"):
                task = t.attrib
                for k, v in task.items():
                    task[k] = get_int(v)
                rts.append(task)

        return rts_id, rts

    #context = et.iterparse(xmlfile, events=('start','end',))
    context = et.iterparse(xml_file, events=('start',))
    context = iter(context)
    event, root = context.__next__()
    #event, root = context.next()

    results = []

    for event, elem in context:
        rts_id, rts = get_rts_from_element(elem)

        if rts:
            # check fu
            fus = [task["C"] / task["T"] for task in rts]
            fu = sum(fus)
            max_percent = (fu * max_fu) / 100
            for task_fu in fus:
                if task_fu > max_percent:
                    break;

            # evaluate the task set
            ok_sched, ok_wcrt, result = test_methods(rts)

            if not ok_sched:
                print("Test results are not equal for task-set {0} (sched).".format(rts_id))
            if not ok_wcrt:
                print("Test results are not equal for task-set {0} (wcrt).".format(rts_id))

            results.append((rts_id, result))

        root.clear()
    del context

    return results


def process_results(results, file):
    sched_results = {}
    sched_results["rta"] = [0, 0, 0]  # sched, non-sched
    sched_results["rta-d"] = [0, 0, 0]  # sched, non-sched
    sched_results["rta-di"] = [0, 0, 0]  # sched, non-sched
    sched_results["rta-het"] = [0, 0, 0]  # sched, non-sched
    sched_results["rta-sjodin"] = [0, 0, 0]  # sched, non-sched
    sched_results["rta-3"] = [0, 0, 0]  # sched, non-sched
    sched_results["rta-3alt"] = [0, 0, 0]  # sched, non-sched

    result_list_dict = []

    for rts_id, result in results:
        for k, v in result.items():
            if v[0]:
                sched_results[k][0] += 1
            else:
                sched_results[k][1] += 1
            sched_results[k][2] += sum(v[2])
            for t, ceils in enumerate(v[2]):
                result_list_dict.append({"method": k, "rts": rts_id, "task": t, "ceils": ceils})

    print(sched_results)

    # Generate Pandas DataFrame for easy data analysis
    result_df = pd.DataFrame(result_list_dict)

    import os
    dump_file = os.path.splitext(file)[0] + '_df.dump'
    pdf_file = os.path.splitext(file)[0] + '.pdf'

    # Plot
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig, ax = plt.subplots(1, 1)
    ax.margins(0.5, 0.5)
    ax.set_ylabel("avg. ceils")
    ax.set_xlabel("task")
    labels = []
    for m, mg in result_df.groupby(["method"]):
        labels.append(m)
        mg.groupby(["task"])[["ceils"]].mean().plot(ax=ax)
    ax.legend(labels, numpoints=1, loc="best", prop={'size': 9})
    plt.savefig(pdf_file, bbox_inches="tight")
    plt.close(fig)

    # Save DataFrame to file
    with open(dump_file, "wb") as outfile:
        pickle.dump(result_df, outfile)


def test_methods(rts):
    results = {}

    results["rta"] = rta_joseph(rts)
    results["rta-d"] = rta_davis(rts)
    results["rta-di"] = rta_davis_inc(rts)
    results["rta-het"] = rta_het(rts)
    results["rta-sjodin"] = rta_sjodin(rts)
    results["rta-3"] = rta3(rts)
    results["rta-3alt"] = rta3_alt(rts)

    # check if all the methods give the same result
    sched_result = []
    for k, v in results.items():
        sched_result.append(v[0])
    sched = sched_result.count(sched_result[0]) == len(sched_result)

    # check that all the methods give the same wcrt for each task
    # same_wcrt = True
    # if not sched:
    #     wcrt = results["rta"][1]
    #     for k, v in results.items():
    #         for idx, w in enumerate(wcrt):
    #             if w != v[1][idx]:
    #                 return sched, False, results

    return sched, True, results


def rta3(rts):
    wcrt = [0] * len(rts)
    ceils = [0] * len(rts)
    a = [0] * len(rts)
    i = [0] * len(rts)
    schedulable = True

    for idx, task in enumerate(rts):
        a[idx] = task["C"]
        i[idx] = task["T"]

    t = rts[0]["C"]
    wcrt[0] = rts[0]["C"]

    for idx, task in enumerate(rts[1:], 1):
        t_mas = t + task["C"]

        while schedulable:
            t = t_mas

            for jdx, jtask in enumerate(rts[:idx]):
                if t_mas > i[jdx]:
                    tmp = math.ceil(t_mas / jtask["T"])
                    a_tmp = tmp * jtask["C"]
                    t_mas += (a_tmp - a[jdx])
                    ceils[idx] += 1

                    if t_mas > task["D"]:
                        schedulable = False
                        break

                    a[jdx] = a_tmp
                    i[jdx] = tmp * jtask["T"]

            if t == t_mas:
                break

        wcrt[idx] = t

        if not schedulable:
            break

    return [schedulable, wcrt, ceils]


def rta3_alt(rts):
    """ Alternative implementation """
    wcrt = [0] * len(rts)
    ceils = [0] * len(rts)
    a = [0] * len(rts)
    i = [0] * len(rts)
    schedulable = True

    for idx, task in enumerate(rts):
        a[idx] = task["C"]
        i[idx] = task["T"]

    r_prev = 0

    for idx, task in enumerate(rts):
        r = r_prev + task["C"]

        while (r > r_prev) and (r <= task["D"]):
            r_prev = r

            for jdx, t in enumerate(rts[:idx]):
                if r > i[jdx]:
                    tmp = math.ceil(r / t["T"])
                    a_tmp = tmp * t["C"]
                    r += (a_tmp - a[jdx])
                    ceils[idx] += 1

                    if r > task["D"]:
                        break

                    a[jdx] = a_tmp
                    i[jdx] = tmp * t["T"]

        wcrt[idx] = r

        if r > task["D"]:
            schedulable = False
            break

    return [schedulable, wcrt, ceils]


def rta_sjodin(rts):
    """ calcula el wrct de cada tarea del str por sjodin """
    ceils = [0] * len(rts)
    wcrt = [0] * len(rts)
    schedulable = True
    r = 0

    for i, task in enumerate(rts, 0):
        c, d = task["C"], task["D"]
        r = r + c

        while True:
            w = 0

            for taskp in rts[:i]:
                cp, tp = taskp["C"], taskp["T"]
                w += math.ceil(r / tp) * cp
                ceils[i] += 1

            w = c + w

            if r == w:
                break

            r = w

            if r > d:
                schedulable = False
                break

        wcrt[i] = r

        if not schedulable:
            break

    return [schedulable, wcrt, ceils]


def rta_het(rts):
    """ calcula el wcrt de cada tarea segun el metodo het """
    def workload(i, b, n):
        if i < 0:
            return 0

        ci = rts[i]["C"]
        ti = rts[i]["T"]

        if b <= last_psi[i]:
            return last_workload[i]

        f = math.floor(b / ti)
        c = math.ceil(b / ti)
        ceils[n] += 2

        branch0 = b - f * (ti - ci) + workload(i - 1, f * ti, n)
        branch1 = c * ci + workload(i - 1, b, n)

        last_psi[i] = b
        last_workload[i] = min(branch0, branch1)

        return last_workload[i]

    ceils = [0] * len(rts)
    wcrt = [0] * len(rts)
    last_workload = [0] * len(rts)
    last_psi = [0] * len(rts)

    schedulable = True

    for i, task in enumerate(rts, 0):
        c = task["C"]
        d = task["D"]
        w = workload(i - 1, d, i)

        if w + c > d:
            schedulable = False
            break

        wcrt[i] = w + c

    return [schedulable, wcrt, ceils]


def rta_davis_inc(rts):
    ceils = [0] * len(rts)
    wcrt = [0] * len(rts)
    inter = [0] * len(rts)
    schedulable = True

    for idx, task in enumerate(rts):
        r_prev = task["C"]
        r = task["C"]

        for jdx, jtask in enumerate(rts[:idx]):
            inter[jdx] = math.ceil(r_prev / jtask["T"]) * jtask["C"]
            ceils[idx] += 1
            r += inter[jdx]

        while (r > r_prev) and (r <= task["D"]):
            r_prev = r

            for jdx, t in enumerate(rts[:idx]):
                tmp = math.ceil(r / t["T"]) * t["C"]
                r += tmp - inter[jdx]
                inter[jdx] = tmp
                ceils[idx] += 1

        wcrt[idx] = r

        if r > task["D"]:
            schedulable = False
            break

    return [schedulable, wcrt, ceils]


def rta_davis(rts):
    ceils = [0] * len(rts)
    wcrt = [0] * len(rts)
    schedulable = True

    for idx, task in enumerate(rts):
        r_prev = 0
        r = task["C"]

        while (r > r_prev) and (r <= task["D"]):
            r_prev = r
            r = task["C"]

            for t in rts[:idx]:
                r += math.ceil(r_prev / t["T"]) * t["C"]
                ceils[idx] += 1

        wcrt[idx] = r

        if r > task["D"]:
            schedulable = False
            break

    return [schedulable, wcrt, ceils]


def rta_joseph(rts):
    """ Joseph and Pandya RTA analysis """
    wcrt = [0] * len(rts)
    ceils = [0] * len(rts)
    schedulable = True

    wcrt[0] = rts[0]["C"]

    for i, task in enumerate(rts[1:], 1):
        r = 1
        c, t, d = task["C"], task["T"], task["D"]

        while schedulable:
            w = 0

            for task_p in rts[:i]:
                cp, tp = task_p["C"], task_p["T"]
                w += math.ceil(r / tp) * cp
                ceils[i] += 1

            w = c + w

            if r == w:
                break

            r = w

            if r > d:
                schedulable = False

        wcrt[i] = r

        if not schedulable:
            break

    return [schedulable, wcrt, ceils]


def get_args():
    """ Command line arguments """
    parser = ArgumentParser(description="Evaluate schedulability tests.")
    parser.add_argument("file", help="XML file with RTS", type=str)
    parser.add_argument("fu", help="FU upper bound per task", type=int, default=30)
    return parser.parse_args()


def main():
    args = get_args()

    if not os.path.isfile(args.file):
        print("{0}: file not found.".format(args.file))
        sys.exit(1)

    if args.fu <= 0 or args.fu > 100:
        print("Invalid fu: {0}.".format(args.fu))
        sys.exit(1)

    # fu and number of tasks from the filename
    rts_data = [int(s) for s in re.findall('\d+', ntpath.basename(args.file))]

    result = test_rts(args.file, args.fu)
    process_results(result, args.file)


if __name__ == '__main__':
    main()
