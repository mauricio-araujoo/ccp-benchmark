
#ifndef EVENTS_H
#define EVENTS_H

#include <unordered_map>
#include <memory>

namespace ccp {

template <class EventType, class BaseActionType>
class Events {
public:

    template <class ActionType> requires std::is_base_of_v<BaseActionType, ActionType>
    std::weak_ptr<ActionType> add_action(EventType event) {
        auto ptr = std::make_shared<ActionType>();
        actions_[event].push_back(ptr);
        return ptr;
    }

    template <class ActionType> requires std::is_base_of_v<BaseActionType, ActionType>
    std::weak_ptr<ActionType> add_action(EventType event, ActionType&& action) {
        auto ptr = std::make_shared<ActionType>(std::move(action));
        actions_[event].push_back(ptr);
        return ptr;
    }

    template <typename... Args>
    void notify(EventType event, Args... args) {
        for (auto& action : actions_[event]) {
            action->notify(args...);
        }
    }

    void clear() {
        actions_.clear();
    }

private:
    std::unordered_map<EventType, std::vector<std::shared_ptr<BaseActionType>>> actions_;

};

} // ccp

#endif //EVENTS_H
